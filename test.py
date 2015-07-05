from libmihmodule import MihasherPython
t = MihasherPython( 10, 20)
t.setK(23)
t.load_bin_codes('data/lsh_64_sift_1M.mat', 'B')

