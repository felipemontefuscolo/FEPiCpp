// This file is part of FEPiC++, a toolbox for finite element codes.
//
// FEPiC++ is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// FEPiC++ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// FEPiC++. If not, see <http://www.gnu.org/licenses/>.

#ifndef FEPIC_IMESHNAMEHANDLER_HPP
#define FEPIC_IMESHNAMEHANDLER_HPP

class _iMeshNameHandler
{
protected:
  
  _iMeshNameHandler() : _is_family(false) {}
  
public:  
  /** @param a file name (without extension) or a path in the form <tt>foo/bar/</tt>
   *  @example if <tt>foo</tt> and <tt>bar</tt> are two directories:
   *           -foo/file  (ok)
   *           -file.vtk  (bad)
   *           -foo/bar   (error)
   *           -foo/bar/  (ok)
   *           -foo/bar/jow (ok only if <tt>bar</tt> exists)
   */ 
  void setOutputFileName(std::string const& name)
  {
    _out_basename = getBaseName(name);
    _out_path     = getRelativePath(name);
    _out_extension= getExtension(name);
    
    if (_out_basename.empty())
      _out_basename = "untitled";
    
    struct stat  dir_stat;
   
    /* check if path passed by user exists */
    bool exists = static_cast<bool>(!lstat(_out_path.data(), &dir_stat));
    bool isdir  = S_ISDIR(dir_stat.st_mode);
    
    if (exists && isdir)
      return;
    else if (!exists)
    {
      // try to create
      if (mkdir(_out_path.data(), S_IRWXU | S_IRGRP | S_IXGRP | S_IXOTH))
      {
        std::cout << "ERROR: can not create directory `"<<_out_path<<"'"<< std::endl;
        throw;
      }
      return;
    }
    else // so _out_path_ exists but it's not a path ... strange, ahn?
    {
      std::cout << "ERROR: can not guess path name, sorry ..." << std::endl;
      throw;
    }
    
  }

  
protected:  
  /// TODO
  void printNames()
  {
    std::cout << "Input file: "               << std::endl;
    std::cout << "meshfile: " << _in_meshfile << std::endl;
    std::cout << "basename: " << _in_basename << std::endl;
    std::cout << "suffix:   " << _in_extension<< std::endl;
    std::cout << "path:     " << _in_path     << std::endl;
    std::cout << "Output file: "              << std::endl;
    
  }

protected:

  /* 
   * Before read a mesh, this function must be called
   * @param filename file name
   * @param extension expected extension
   * @param is_family output
   */ 
  bool _registerFile(std::string filename, std::string const& extension)
  {
    std::ifstream file(filename);
    
    FEPIC_ASSERT(!file.fail(), "can not find mesh file", std::invalid_argument);
    file.close();
    
    _in_meshfile  = ::stripTrailingSpaces(filename);
    _in_extension = ::getExtension(filename);
    _in_basename  = ::getBaseName(filename);
    _in_path      = ::getRelativePath(filename);
    
    if (_in_extension.empty())
      std::cout << "WARNING: mesh file without extension" << std::endl;
    else if ( _in_extension != extension)
    {
      std::cout << "WARNING: wrong file extension: "
                << extension << " expected, got " << _in_extension << std::endl;
    }
    
    /* default values */
    _out_path     = _in_path;
    _out_basename = _in_basename; 
    
    return 0;
  }
 
  std::string _popNextName(int filenum, std::string const& ext)
  {
    // filenum : a suffix to basename; the series number
    
    if (_is_family)
      return _out_path+_out_basename+itoafill0(filenum, FEPIC_FILE_FILL)+ext;
    else
      return _out_path+_out_basename+ext;
  }

  std::string _in_meshfile;  // eg.   /home/user/test.msh
  std::string _in_basename;  // eg.   test
  std::string _in_extension; // eg.   .msh
  std::string _in_path;      // eg.   /home/user/
  std::string _out_basename; // eg.   result_file
  std::string _out_path;     // eg.   /home/user/result/
  std::string _out_extension;// eg.
  bool        _is_family; // se o output sera impresso como familia  

};



#endif
