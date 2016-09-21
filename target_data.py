import numpy as np
import config


class TargetData:
    # target_data provides a list of targets used to specify RA and Dec for a given .fits file's target object.
    # NOTE: data in this array is of type Byte.
    target_data = np.loadtxt(fname=config.target_data + 'target_data.txt',
                             dtype=bytes,
                             delimiter=',')

    def __init__(self):
        """
        Initiates TargetData class
        """
        pass

    def bytes_to_str(self, list_input):
        """
        Converts data type of imported array elements from
        bytes to string.

        :type list_input: list
        """
        for i in range(len(list_input)):
            x = list_input[i].decode('utf-8')
            list_input[i] = x
        return list_input

    def target_dict(self, name_list, val_list):
        """
        Takes list of names and list of values and returns a
        dictionary

        :type name_list: list used for dictionary key
        :type val_list: list used for dictionary value for each key
        """
        # lower all strings in name_list for consistency with comparisons in photometry.net
        for i in range(len(name_list)):
            name_list[i] = name_list[i].lower()
        result = dict(zip(name_list, val_list))
        return result

    def coord_lookup(self, file1, dict1):
        """
        Used in script_loop() to relate parsed file name to a given coordinate. Filename is split according to
        conditionals below.

        :type file1: string representing filename that will be compared to keys of dict1
        :type dict1: dictionary that is used to lookup coordinates via the keyvalue to file1

        """
        i = ''
        file1_rmus = file1.replace('_', '', 1)
        file1_nows = file1_rmus.replace(' ', '', 1)
        for key in dict1.keys():
            if key.lower() in file1.lower():
                i = key
            elif key.lower() in file1_rmus.lower():
                i = key
            elif key.lower() in file1_nows.lower():
                i = key
            else:
                pass
        coord = dict1.get(i)
        return coord

    """
    # All values from target_data are Byte data type.
    target_name = target_data[:, 0].tolist()
    target_RA = target_data[:, 1].tolist()
    target_dec = target_data[:, 2].tolist()
    # Converts Byte data in each list to string data type.
    target_name = bytes_to_str(target_name)
    target_RA = bytes_to_str(target_RA)
    target_dec = bytes_to_str(target_dec)
    # create dictionaries of RA and dec for each respective target
    RA_dict = target_dict(target_name, target_RA)
    dec_dict = target_dict(target_name, target_dec)
    """