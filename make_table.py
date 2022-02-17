import pandas as pd
from sys import argv

file = argv[1]

df = pd.read_csv(file,skipinitialspace=True)

if len(argv)==3:
    print('# %s  '%argv[2])

print('| Sensor 1 | Sensor 2 | Profiles | n | r2 | a (line) | b (line, 10<sup>-4</sup> m<sup>-1</sup>) | b (offset, 10<sup>-4</sup> m<sup>-1</sup>) |  ')
print('|:---------|:---------|:---------|:--|:---|:--|:--|:--|  ')
for index,row in df.iterrows():
    print('| %s | %s | %d | %d | %.2f | %.2f &pm; %.2f | %.1f &pm; %.1f | %.1f &pm; %.1f |  '%
          (row['inst1'],row['inst2'],row['num_profiles'],row['n'],row['rsq'],row['line_a'],
           row['line_a_95CI'],row['line_b']*1e4,row['line_b_95CI']*1e4,row['offset_b']*1e4,row['offset_b_95CI']*1e4))
