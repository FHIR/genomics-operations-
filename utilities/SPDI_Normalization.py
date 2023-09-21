from glob import glob
from math import ceil
from pathlib import Path
from threading import Lock

MAX_REFSEQ_CACHE_SIZE = 300 * 1024 * 1024  # 450MB


def hex2Code(code):
    """
    Convert a 4 byte integer to a sequence code. See `code2Hex()` for details
    """
    match code:
        case 0x0: return 'A'
        case 0x1: return 'C'
        case 0x2: return 'G'
        case 0x3: return 'T'
        case 0x4: return 'U'
        case 0x5: return 'R'
        case 0x6: return 'Y'
        case 0x7: return 'K'
        case 0x8: return 'M'
        case 0x9: return 'S'
        case 0xA: return 'W'
        case 0xB: return 'B'
        case 0xC: return 'D'
        case 0xD: return 'H'
        case 0xE: return 'V'
        case 0xF: return 'N'
        case _:
            raise NotImplementedError(f"unsupported code '{code}'")


class RefSeq:
    def __init__(self, refSeq, length):
        """
        Store the length so we can easily tell when the index operator receives an out of bounds subscript because we
        want to represent each symbol using 4 bits and we have 16 symbols, so there's no room left to represent an
        invalid symbol which denotes end-of-sequence.
        """
        self.__packedRefSeq = refSeq
        self.__length = length

        # Size in bytes
        self.size = ceil(length / 2)

    def __getitem__(self, subscript):
        if isinstance(subscript, slice):
            codes = []
            for i in range(*subscript.indices(self.__length)):
                codes.append(self[i])
            return "".join(codes)

        byte = self.__packedRefSeq[subscript // 2]

        if subscript % 2 == 0:
            return hex2Code((byte >> 4) & 0xF)

        return hex2Code(byte & 0xF)


class RefSeqCache:
    def __init__(self):
        self.cache = {}
        self.size = 0

    def __getKey(self, build, refSeqName):
        return f'{build}/{refSeqName}'

    def contains(self, build, refSeqName):
        return self.__getKey(build, refSeqName) in self.cache

    def get(self, build, refSeqName):
        return self.cache[self.__getKey(build, refSeqName)]

    def add(self, build, refSeqName, refSeq):
        """
        Evacuates enough items from the cache until we have room for the new one and then adds it to the cache.
        """

        while self.size + refSeq.size > MAX_REFSEQ_CACHE_SIZE:
            key, cachedRefSeq = next(iter(self.cache.items()))
            self.size -= cachedRefSeq.size
            del self.cache[key]

        self.cache[self.__getKey(build, refSeqName)] = refSeq
        self.size += refSeq.size


refSeqCache = RefSeqCache()
refSeqCacheLock = Lock()


def get_ref_seq(build, refSeqName):
    if not refSeqCache.contains(build, refSeqName):
        try:
            refSeqFilePattern = f'refseq/{build}/{refSeqName}_*.seq'
            foundFiles = glob(refSeqFilePattern)
            if len(foundFiles) != 1:
                raise ValueError(f'failed to find expected refseq file for build "{build}" and RefSeq "{refSeqName}" ')

            refSeqFile = foundFiles[0]
            refSeqLength = Path(refSeqFile).stem.rpartition('_')[-1]
            with open(refSeqFile, mode='rb') as file:
                refSeq = RefSeq(file.read(), int(refSeqLength))
        except Exception as err:
            print(f"failed to read refseq file: {err=}, {type(err)=}")
            raise
        refSeqCache.add(build, refSeqName, refSeq)
    return refSeqCache.get(build, refSeqName)


def get_normalized_spdi(ref_seq, pos, ref, alt, build):
    # Need to serialise this if we can't keep all the ref seq data in memory
    with refSeqCacheLock:
        return get_normalized_spdi_impl(ref_seq, pos, ref, alt, build)


def get_normalized_spdi_impl(ref_seq, pos, ref, alt, build):
    ref_seq_fasta = get_ref_seq(build, ref_seq)

    # Step 0
    start = pos
    end = pos + len(ref)

    # Step 1
    # a. Trim Common Suffix

    small = len(alt) if len(alt) <= len(ref) else len(ref)

    suffix_length = 0

    last = -1
    while small != 0:
        if ref[last] == alt[last]:
            suffix_length += 1
            end -= 1
            last -= 1
            small -= 1
        else:
            break

    if suffix_length > 0:
        ref = ref[0:-suffix_length]
        alt = alt[0:-suffix_length]

        # b. Trim Common Prefix

    small = len(alt) if len(alt) <= len(ref) else len(ref)

    prefix_length = 0

    first = 0
    while small != 0:
        if ref[first] == alt[first]:
            prefix_length += 1
            start += 1
            first += 1
            small -= 1
        else:
            break

    if prefix_length > 0:
        ref = ref[prefix_length:]
        alt = alt[prefix_length:]

    # Step 2
        # a. Check if both are empty
        # Not Possible in this case

        # b. Check if both are non-empty

    if ref and alt:
        return f"{ref_seq}:{start}:{ref}:{alt}"

        # c. Check if one is empty

    if alt:
        allele = alt
        type_ = 'I'
    else:
        allele = ref
        type_ = 'D'

    # Step 3
        # a. Left Roll

    if type_ == 'I':
        left_roll_bound = start

        temp_allele = allele

        while (str(ref_seq_fasta[left_roll_bound-1]).upper() == temp_allele[-1].upper()):
            temp_allele = str(ref_seq_fasta[left_roll_bound-1]).upper() + temp_allele[:-1]
            left_roll_bound -= 1

            # . b. Right Roll

        right_roll_bound = start

        temp_allele = allele

        while (str(ref_seq_fasta[right_roll_bound]).upper() == temp_allele[0].upper()):
            temp_allele = temp_allele[1:] + str(ref_seq_fasta[right_roll_bound]).upper()
            right_roll_bound += 1

        # Step 4
            # a. Prepend
        if ref_seq_fasta[left_roll_bound:start]:
            ref = str(ref_seq_fasta[left_roll_bound:start]) + ref
            alt = str(ref_seq_fasta[left_roll_bound:start]) + alt

            # b. append
        if ref_seq_fasta[start:right_roll_bound]:
            ref = ref + str(ref_seq_fasta[start:right_roll_bound])
            alt = alt + str(ref_seq_fasta[start:right_roll_bound])

            # c. Modify Start and End, and Return SPDI
        start = left_roll_bound
        end = right_roll_bound

        return (f"{ref_seq}:{start}:{ref.upper()}:{alt.upper()}")

    else:
        left_roll_bound = start

        while (str(ref_seq_fasta[left_roll_bound-1]).upper() == str(ref_seq_fasta[left_roll_bound + (len(allele)-1)]).upper()):
            left_roll_bound -= 1

            # . b. Right Roll

        right_roll_bound = start

        while (str(ref_seq_fasta[right_roll_bound]).upper() == (ref_seq_fasta[right_roll_bound + (len(allele))]).upper()):
            right_roll_bound += 1

        # Step 4
            # a. Prepend
        if ref_seq_fasta[left_roll_bound:start]:
            ref = str(ref_seq_fasta[left_roll_bound:start]) + ref
            alt = str(ref_seq_fasta[left_roll_bound:start]) + alt

            # b. append
        if ref_seq_fasta[start:right_roll_bound]:
            ref = ref + str(ref_seq_fasta[start:right_roll_bound])
            alt = alt + str(ref_seq_fasta[start:right_roll_bound])

            # c. Modify Start and End, and Return SPDI
        start = left_roll_bound
        end = right_roll_bound

        return (f"{ref_seq}:{start}:{ref.upper()}:{alt.upper()}")
