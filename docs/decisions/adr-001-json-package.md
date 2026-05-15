# json package

## Context and Problem Statement

Python provides several different packages to parse json files that differ in performance:

package | documentation
---|---
json | https://docs.python.org/3/library/json.html
ijson | https://pypi.org/project/ijson/
msgspec | https://pypi.org/project/msgspec/
orjson | https://pypi.org/project/orjson/
ujson | https://pypi.org/project/ujson/

`json` is the standard package. We would like to know how the other packages perform in comparison when reading a big json file. After installing all packages we were using the following script to test performance:

```python
import ijson
import json
import msgspec
import orjson
import time
import ujson

file_path = 'path/to/big.json'

# json
start_time = time.perf_counter()
with open(file_path, 'r') as file:
    data = json.load(file)
    print(data['header'])
end_time = time.perf_counter()
print(f"json: {end_time - start_time}")

# ijson
start_time = time.perf_counter()
with open(file_path, 'r') as file:
    for item in ijson.items(file, 'header'):
        print(item)
end_time = time.perf_counter()
print(f"ijson: {end_time - start_time}")

# msgspec
start_time = time.perf_counter()
with open(file_path, 'r') as file:
    data = msgspec.json.decode(file.read())
    print(data['header'])
end_time = time.perf_counter()
print(f"msgspec: {end_time - start_time}")

# orjson
start_time = time.perf_counter()
with open(file_path, 'r') as file:
    data = orjson.loads(file.read())
    print(data['header'])
end_time = time.perf_counter()
print(f"orjson: {end_time - start_time}")

# ujson
start_time = time.perf_counter()
with open(file_path, 'r') as file:
    data = ujson.loads(file.read())
    print(data['header'])
end_time = time.perf_counter()
print(f"ujson: {end_time - start_time}")
```

Both `msgspec` and `orjson` perform the fastest.

## Decision

## Consequences

### Positive

-

### Negative

-

### Neutral
