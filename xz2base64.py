import base64

"""
To insert the base64 file in the code:

import base64, lzma

b64_file = """ """

compressed_b64_file = base64.b64decode(b64_file)
file = lzma.decompress(compressed_b64_file)

"""


with open("compressed_file.xz", "rb") as f:
    encoded = base64.b64encode(f.read()).decode()

with open("compressed_file.base64.py", "w") as f:
    f.write(f'eggNOG-index = """{encoded}"""')

