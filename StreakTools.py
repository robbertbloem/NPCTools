from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

import numpy
import matplotlib 
import matplotlib.pyplot as plt

import struct

def import_streak(paf):

    with open(paf, 'rb') as f:
        data = f.read()

    comment_length = struct.unpack_from("h", data, offset = 2)[0]
    width = struct.unpack_from("h", data, offset = 4)[0]
    height = struct.unpack_from("h", data, offset = 6)[0]
    x_offset = struct.unpack_from("h", data, offset = 8)[0]
    y_offset = struct.unpack_from("h", data, offset = 10)[0]
    type = struct.unpack_from("h", data, offset = 12)[0]

#     print(comment_length, width, height, x_offset, y_offset, type)

    comment_string = str(data[64:comment_length+64])
#     comment_string = comment_string.replace(r"\r\n", ",")
#     comment = comment_string.split(",")
#     for c in comment:
#         print(c)

    data2 = struct.unpack_from("i" * width * height, data, offset = comment_length+64)

    data2 = numpy.reshape(data2, (height, width))

    x_axis_offset = (comment_length + 64) + width * height * 4
    y_axis_offset = x_axis_offset + width * 4

    x_axis = struct.unpack_from("f" * width, data, offset = x_axis_offset) 
    y_axis = struct.unpack_from("f" * height, data, offset = y_axis_offset)

    x_axis = numpy.array(x_axis)
    y_axis = numpy.array(y_axis)

    return data2, x_axis, y_axis, comment_string

if __name__ == "__main__": 
    pass