import h5py

def list_datasets(file_path):
    '''printing all datasets in an HDF5 file'''
    with h5py.File(file_path, 'r') as f:
        def print_name(name, obj):
            if isinstance(obj, h5py.Dataset):
                print(name)
        f.visititems(print_name)

file_path = '/geos/netdata/oosa/assignment/lvis/2015/ILVIS1B_AQ2015_1012_R1605_067095.h5'
list_datasets(file_path)
