import tests.utilities as tu
import json


"""
Test Json Sorting Tests
-------------------------------------
"""
simple = {
    "type": "simple-fruit",
    "example": "apple",
    "seeds?": "yes"
}

simple_nest = {
   "fruits": [
      {
         "type": "simple-fruits",
         "examples": [
            {
               "name": "apple"
            },
            {
               "name": "cherries"
            }
         ]
      },
      {
         "type": "aggregate",
         "examples": [
            {
               "name": "drupelets",
               "members": [
                {
                   "common-name": "raspberry"
                },
                {
                    "common-name": "blackberry"
                }
               ]
            },
            {
               "name": "achenes",
               "common-name": "strawberry"
            }
         ]
      }
   ]
}


def test_sorting_json_simple():
    sorted_actual = tu.sort_output_jsons(simple)
    expected_loc = f'{tu.JSON_SORT_OUTPUT_DIR}simple.json'
    with open(expected_loc) as expected_file:
        expected_json = json.load(expected_file)
        assert json.dumps(sorted_actual) == json.dumps(expected_json)


def test_sorting_json_nested():
    sorted_actual = tu.sort_output_jsons(simple_nest)
    expected_loc = f'{tu.JSON_SORT_OUTPUT_DIR}nested.json'
    with open(expected_loc) as expected_file:
        expected_json = json.load(expected_file)
        assert json.dumps(sorted_actual) == json.dumps(expected_json)
