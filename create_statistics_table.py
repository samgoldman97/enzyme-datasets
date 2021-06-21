"""Create dataset statistics table"""

import os
import pandas as pd

def main(): 
    processed_data_dir = "data/processed" 

    file_names = os.listdir(processed_data_dir)
    file_names = [i for i in file_names if i.endswith(".csv")]
    dataset_summary_list = []
    for file_name in file_names: 
        df = pd.read_csv(os.path.join(processed_data_dir, file_name),  
                         index_col = 0 )

        uniq_seqs = len(pd.unique(df['SEQ']))
        uniq_subs = len(pd.unique(df['SUBSTRATES']))

        total_pairs = len(df.groupby(["SEQ", "SUBSTRATES"]))
        activity_col = set(df.columns.values.tolist()).difference(["SEQ", "SUBSTRATES"])
        activity_col = list(activity_col)[0]
        uniq_vals = len(pd.unique(df[activity_col]))
        if uniq_vals == 2: 
            dataset_type = "Binary"
        elif uniq_vals == 3: 
            dataset_type = "Categorical"
        elif uniq_vals > 3: 
            dataset_type = "Regression"
        else:
            raise ValueError()

        new_entry = {"Dataset Name" : file_name,
                     "Dataset Type" : dataset_type, 
                     "Dataset Class" : file_name.split("_")[0].split(".")[0], 
                     "# Enz." : uniq_seqs, 
                     "# Sub." : uniq_subs, 
                     "# Pairs" : total_pairs}
        dataset_summary_list.append(new_entry)

    df = pd.DataFrame(dataset_summary_list)
    df = df.reset_index(drop=True)
    df = df.set_index("Dataset Class")
    df.sort_values(by="# Pairs", inplace=True)

    markdown_str = df.to_markdown()
    open("markdown_temp.md", "w").write(markdown_str)


if __name__=="__main__": 
    main()
