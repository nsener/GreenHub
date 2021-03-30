import json
import pandas as pd
# load json file
with open("det_carb_lc.json", 'r' ) as f:
    data = json.load(f)

# convert json file content(list of dictinaries) to pandas dataframe
df = pd.DataFrame(data)
# get keys of each opened hub location and convert them to a list
# and append the result inside very same dataframe under 'OpenedHubIds' column
df['OpenedHubIds'] = df['OpenedHubLocs'].apply(lambda x: list(x.keys()))
print(df['OpenedHubIds'])