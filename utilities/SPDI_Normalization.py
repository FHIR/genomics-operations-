from glob import glob
from math import floor
from pathlib import Path
from threading import Lock


def hex_to_code(hex):
    """
    Convert a 3 bit integer to a sequence code. See `code_to_hex()` in pack_refseq.py for details.
    """
    match hex:
        case 0x0: return 'A'
        case 0x1: return 'C'
        case 0x2: return 'G'
        case 0x3: return 'T'
        case 0x4: return 'N'
        # Blow up if we get any unexpected hex-encoded code
        case _:
            raise NotImplementedError(f"unexpected hex-encoded code: '{hex}'")


class RefSeq:
    # TODO: Consider adding `__str__()` method (and maybe `__repr__()` too?)

    def __init__(self, packed_ref_seq, len):
        self.__packed_ref_seq = packed_ref_seq
        # Store the RefSeq length so we can easily tell when the index operator receives an out of bounds subscript
        # because we want to represent each code using 3 bits and we don't want to introduce an end-of-sequence code, in
        # case we'll want to leverage the unused codes for something else.
        self.__len = len

    def __getitem__(self, subscript):
        if isinstance(subscript, slice):
            codes = []
            # TODO: It should be more efficient to process 3 bytes at a time for long ranges.
            for i in range(*subscript.indices(self.__len)):
                codes.append(self[i])
            return "".join(codes)

        if subscript >= self.__len:
            raise IndexError(f"list index out of range: '{subscript}' must be less than '{self.__len}'")

        # Select the byte(s) and the offset of the code we wish to extract from the packed RefSeq data.
        packed_subscript = floor(subscript * 3 / 8)
        data = self.__packed_ref_seq[packed_subscript]
        match subscript % 8:
            case 0: offset = 0
            case 1: offset = 3
            case 2:
                offset = 6
                # This code is packed across two bytes.
                data += self.__packed_ref_seq[packed_subscript + 1] << 8
            case 3: offset = 1
            case 4: offset = 4
            case 5:
                offset = 7
                # This code is packed across two bytes.
                data += self.__packed_ref_seq[packed_subscript + 1] << 8
            case 6: offset = 2
            case 7: offset = 5

        # Move the relevant three bits to the bottom, extract them using a bit mask and convert them back to a code.
        return hex_to_code((data >> offset) & 0x3)


def get_ref_seq(build, ref_seq_name):
    try:
        ref_seq_file_pattern = f'refseq/{build}/{ref_seq_name}_*.refseq'
        found_files = glob(ref_seq_file_pattern)
        if len(found_files) != 1:
            raise ValueError(f'failed to find expected refseq file for build "{build}" and RefSeq "{ref_seq_name}" ')

        ref_seq_file = found_files[0]
        ref_seq_length = int(Path(ref_seq_file).stem.rpartition('_')[-1])

        with open(ref_seq_file, mode='rb') as file:
            refSeq = RefSeq(file.read(), ref_seq_length)
    except Exception as err:
        print(f"failed to read refseq file: {err=}, {type(err)=}")
        raise

    return refSeq


ref_seq_lock = Lock()


def get_ref_seq_subseq(build, ref_seq, start, end):
    # Need to serialise this if we can't keep all the RefSeq data in memory
    with ref_seq_lock:
        return get_ref_seq(build, ref_seq)[start:end]


def get_normalized_spdi(ref_seq, pos, ref, alt, build):
    # Need to serialise this if we can't keep all the RefSeq data in memory
    with ref_seq_lock:
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
