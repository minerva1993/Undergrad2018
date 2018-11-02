import os

path_to_prod = '../Delphes2Flat/output'
print("Looking for files in %s"%path_to_prod)

input_list_file_name = 'file_list.txt'
string_for_processing = ''

#This part is for reconstruction
for file_name in os.listdir(path_to_prod):
  if file_name == ".gitkeep": continue
  tmp_string = ''
  tmp_string += os.path.join(os.path.abspath(path_to_prod), file_name)
  dataset_name = file_name
  tmp_string += ' ./output/hist_' + dataset_name.replace("_sample","")
  string_for_processing += tmp_string + '\n'

with open(input_list_file_name, 'w') as f:
  f.write(string_for_processing)


print("Done")

