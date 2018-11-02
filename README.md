# Undergrad2018
From Delphes2Flat
```{.Bash}
ssh compute-0-2
cd ./path_to/Delphes2Flat
python create_input_file_list.py
cat file_list.txt | xargs -i -P$(nproc) -n2 python run.py
```

Run analyzer
```{.Bash}
root -l run.C'("input.root","outname.root")'
#or
source job.sh
```
