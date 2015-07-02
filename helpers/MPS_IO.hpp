#ifndef _MPS_IO
#define _MPS_IO

#include <sstream>

/**< amend armadillo's hdf5 saving abilities to save several objects into a single HDF5 file */
namespace arma
{

template<typename VT>
bool save_hdf5_binary(const Mat<VT>& mat, const std::string& varname, hid_t file)
{
    bool write_okay = false;
// We need to create a dataset, datatype, and dataspace
    hsize_t dims[2];
    dims[1] = mat.n_rows;
    dims[0] = mat.n_cols;

    hid_t dataspace = arma_H5Screate_simple(2, dims, NULL);   // treat the matrix as a 2d array dataspace
    hid_t datatype  = hdf5_misc::get_hdf5_type<VT>();

    // If this returned something invalid, well, it's time to crash.
    arma_check(datatype == -1, "Mat::save(): unknown datatype for HDF5");

    // MATLAB forces the users to specify a name at save time for HDF5; Octave
    // will use the default of 'dataset' unless otherwise specified, so we will
    // use that.
    hid_t dataset = arma_H5Dcreate(file, varname.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // H5Dwrite does not make a distinction between row-major and column-major;
    // it just writes the memory.  MATLAB and Octave store HDF5 matrices as
    // column-major, though, so we can save ours like that too and not need to
    // transpose.
    herr_t status = arma_H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, mat.mem);
    write_okay = (status >= 0);

    arma_H5Dclose(dataset);
    arma_H5Tclose(datatype);
    arma_H5Sclose(dataspace);
    return write_okay;
}
}

#ifndef SYMS/**< implementation for DenseMPS */
template<typename KT,typename VT>
void save(const MPS);

#else/**< implementation for BlockMPS */
//#include "MPSBlockMat.hpp"

template<typename KT,typename VT>
bool save(const MPSBlockMat<KT,VT>& MPS, const std::string& name)
{
    bool save_okay = false;

    const std::string tmp_name = diskio::gen_tmp_name(name);

    // Set up the file according to HDF5's preferences
    hid_t file = arma_H5Fcreate(tmp_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    for (size_t i=0;i<MPS.size();++i)
    {
        cout<<"saving "<<i<<endl;
        for (const auto& it:MPS[i])
        {
            stringstream tmpname;
//            tmpname << i <<
        }
    }

//    herr_t close_okay = arma_H5Fclose(file);
    save_okay = !(arma_H5Fclose(file)<0);
//    cout<<close_okay<<endl;

    if(save_okay == true) { save_okay = diskio::safe_rename(tmp_name, name); }
    return save_okay;
}
#endif // _USE_SYMMETRIES



#endif // _MPS_IO
