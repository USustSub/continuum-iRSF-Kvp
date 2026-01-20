function save_mat ( file_ , ff , inc ) 

    filename = sprintf(['out/' (file_) '/data/jj%02d.mat'], inc);
    save(filename, 'ff');

end