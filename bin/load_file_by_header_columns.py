import os
import numpy as np

EXTENSION_TO_DELIMITER_MAP = {
    '.txt': '\t',
    '.csv': ',',
}

def get_file_delimiter(file_name):
    '''
        determine which delimiter to use based on the extension
        return the corresponding value from EXTENSION_TO_DELIMITER_MAP
    '''
    ext = os.path.splitext(file_name)[1]

    if not ext in EXTENSION_TO_DELIMITER_MAP:
        raise Exception("Couldn't find extention {0} in map {1}".format(ext, EXTENSION_TO_DELIMITER_MAP))

    return EXTENSION_TO_DELIMITER_MAP[ext]

def homogenize_column_name(string):
    return string.replace(' ', '').replace('_','').upper()

def map_required_columns(headers, required_column_names):
    '''
        Create a map of required header names to their actual indices
        headers : list of string headers that were in the first line
        required_column_names : list of headers that should be in the header list
        raise an exception if any required header is not found
        return a list of indices that correspond to each required header value
    '''
    # make a dict of the required names
    actual_headers = {name.strip().upper(): index for index,name in enumerate(headers)}

    required_indices = []
    for required_header in required_column_names:
        stripped_req_header = homogenize_column_name(required_header)
        if stripped_req_header not in actual_headers:
            raise Exception("Required header \'{0}\' was not found in actual headers {1}".format(required_header, headers))
        required_indices.append(actual_headers[stripped_req_header])
    return required_indices

def load_file_by_headers(file_name, required_column_names, required_column_converters):
    '''
        Load a file and return data for columns specified by their corresponding headers.
        The headers can be specified in any case and spaces or underscores are ignored,
        allowing to pass for example protein_name and match to Protein Name, or ProteinName, etc.
        file_name : path of file to open
        required_column_names : list of strings for columns that should be present in the header
        required_column_converters : list of functions for how to load each desired column
        return a list of dict for each row in the file, where each key in the dict is the required column name,
        and the value is the value from that row.
    '''
    if not isinstance(required_column_names, list):
        raise Exception('A list of column names must be given')
    if not isinstance(required_column_converters, list):
        raise Exception('A list of column converters must be given')
    if len(required_column_names) != len(required_column_converters):
        raise Exception('The required names and converter lists must be the same size')

    with open(file_name, 'r') as f:
        # get the delimiter for this file type
        delimiter = get_file_delimiter(file_name)
        # get the header line and split it
        headers = [homogenize_column_name(x) for x in f.readline().split(delimiter)]
        # get the map of columns we want to their index
        req_header_indices = map_required_columns(headers, required_column_names)
        
        # create the map of column index to converter function
        converters = {index: func for index,func in zip(req_header_indices, required_column_converters)}
        # load the file 
        data = np.genfromtxt(file_name, delimiter=delimiter, skip_header=1, dtype=None, encoding=None, converters=converters, usecols=req_header_indices)
        
        output_items = []
        num_cols = len(required_column_names)
        if num_cols == 1:
            # if there's only one column, then we have to parse the output a little differently
            key = required_column_names[0]
            for value in data:
                output_items.append({key: value})
        else:
            for row in data:
                item = {}
                for name,value in zip(required_column_names, row):
                    item[name] = value
                output_items.append(item)
    return output_items