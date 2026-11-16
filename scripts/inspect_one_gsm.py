#load library
import pandas as pd

#create dataframe for sample GSM to ensure it looks correct
df = pd.read_csv(
    "data/dev/GSM1545535_10_6_5_11.txt.gz",
    sep="\t",
    header=None
)

#view data frame
print(df.head())
print(df.shape)