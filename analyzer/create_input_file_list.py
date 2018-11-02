import os

path_to_prod = '../Delphes2Flat/output'
print("Looking for files in %s"%path_to_prod)

input_list_file_name = 'file_list.txt'
string_for_processing = ''

merge_file = 'job_merge.sh'
string_for_merging = ''

#This part is for reconstruction
for file_name in os.listdir(path_to_prod):
  if file_name == ".gitkeep": continue
  tmp_string = ''
  tmp_string += os.path.join(os.path.abspath(path_to_prod), file_name)
  tmp_string += ' ./output/hist_' + file_name.replace("_sample","")
  string_for_processing += tmp_string + '\n'

with open(input_list_file_name, 'w') as f:
  f.write(string_for_processing)

dataset_list = []
for file_name in os.listdir(path_to_prod):
  if file_name == ".gitkeep": continue
  dataset = file_name.split("_")[0]
  if dataset not in dataset_list:
    dataset_list.append(dataset)

for tmp_set in dataset_list:  
  string_for_merging += 'hadd merged/hist_' + tmp_set + '.root output/hist_' + tmp_set + '*.root\n'

with open(merge_file, 'w') as f:
  f.write(string_for_merging)

print("Done")

