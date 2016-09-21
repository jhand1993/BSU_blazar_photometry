import os

# creating variables
pipeline_root = '/Users/jaredhand/challis_research/blazar_photometry/'

astrometry = pipeline_root + 'astrometry/'
stacking = pipeline_root + 'stacking/'
sex = pipeline_root + 'sex/'
catalogs = pipeline_root + 'catalogs/'
finished_stacked = pipeline_root + 'finished_stacked/'
target_data = pipeline_root + 'target_data/'
dirlists = [astrometry, stacking, sex, catalogs, finished_stacked, target_data]


def makedirs():
    # making directories
    for i in dirlists:
        try:
            os.mkdir(i)
        except FileExistsError:
            print(i, 'already exists.')
        except:
            raise

# global variables
obj = 'BLLAC'
r = 5 / 60  # 10 / 3600

if __name__ == '__main__':
    makedirs()
