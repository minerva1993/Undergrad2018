import os

path_to_prod = '/data/users/minerva1993/work/2018_UGLeptoQ/output'
print("Looking for files in %s"%path_to_prod)

input_list_file_name = 'file_list.txt'
string_for_processing = ''

#This part is for reconstruction
for dataset_folder in os.listdir(path_to_prod):
  dataset_path = os.path.join(path_to_prod, dataset_folder)
  for file_name in os.listdir(dataset_path):
    tmp_string = ''
    tmp_string += os.path.join(dataset_path, file_name)
    dataset_name = dataset_folder.replace("_",'')
    tmp_string += ' ' + dataset_name + '_'
    string_for_processing += tmp_string + '\n'

with open(input_list_file_name, 'w') as f:
  f.write(string_for_processing)


print("Done")
