import base64
import sys

xz_file = sys.argv[1]
b64_file = sys.argv[2]

with open(xz_file, "rb") as f:
    encoded = base64.b64encode(f.read()).decode('utf-8')

with open(b64_file, "w") as f:
    f.write(f'XZ_DATA = """{encoded}"""')
