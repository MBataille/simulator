import matplotlib
import shutil
import os

mpl_folder = matplotlib.__file__

MPL_ICONS_FOLDER = '/mpl-data/images/'
mpl_folder = mpl_folder[:mpl_folder.rfind('/')] + MPL_ICONS_FOLDER

ICONS_FOLDER = 'icons/'
files = [f for f in os.listdir(ICONS_FOLDER) if os.path.isfile(os.path.join(ICONS_FOLDER, f))]

for file in files:
	# print(f'moving {ICONS_FOLDER + file} to {mpl_folder + file}')
	shutil.copyfile(ICONS_FOLDER + file, mpl_folder + file)