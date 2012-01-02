#include <cstdio>
#include <sys/stat.h>

#include "meshnamehandler.hpp"
#include "../../util/macros.hpp"
#include "../../util/assert.hpp"



/** @param a file name (without extension) or a path in the form <tt>foo/bar/</tt>
 *  @example if <tt>foo</tt> and <tt>bar</tt> are two directories:
 *           -foo/file  (ok)
 *           -file.vtk  (bad)
 *           -foo/bar   (error)
 *           -foo/bar/  (ok)
 *           -foo/bar/jow (ok only if <tt>bar</tt> exists)
 */ 
void _MeshNameHandler::setOutputFileName(const char* name)
{
  _sofn_already_called = true;
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
      printf("ERROR: can not create directory `%s'\n", _out_path.c_str());
      throw;
    }
    return;
  }
  else // so _out_path_ exists but it's not a path ... strange, ahn?
  {
    printf("ERROR: can not guess path name, sorry ...\n");
    throw;
  }
  
}





/* 
 * Before read a mesh, this function must be called
 * @param filename file name
 * @param extension expected extension
 * @param is_family output
 */ 
bool _MeshNameHandler::_registerFile(std::string filename, std::string const& extension)
{
  
  FILE *file_ptr = fopen(filename.c_str(), "r");
  FEPIC_ASSERT(file_ptr!=NULL, "can not find mesh file", std::invalid_argument);
  fclose(file_ptr);
  
  _in_meshfile  = ::stripTrailingSpaces(filename);
  _in_extension = ::getExtension(filename);
  _in_basename  = ::getBaseName(filename);
  _in_path      = ::getRelativePath(filename);
  
  if (_in_extension.empty())
    printf("WARNING: mesh file without extension\n");
  else if ( _in_extension != extension)
  {
    printf("WARNING: wrong file extension: %s expected, got %s", extension.c_str(), _in_extension.c_str());
  }
  
  /* default values */
  _out_path     = _in_path;
  _out_basename = _in_basename; 
  
  return 0;
}




//const char* _MeshNameHandler::_popNextName(int filenum, std::string const& ext)
//{
  //// filenum : a suffix to basename; the series number
  //if (_is_family)
    //return (_out_path+_out_basename+itoafill0(filenum, FEPIC_FILE_FILL)+ext).c_str();
  //else
    //return (_out_path+_out_basename+ext).c_str();
//}

std::string _MeshNameHandler::_popNextName(int filenum, std::string const& ext)
{
  // filenum : a suffix to basename; the series number
  if (_is_family)
    return _out_path+_out_basename+itoafill0(filenum, FEPIC_FILE_FILL)+ext;
  else
    return _out_path+_out_basename+ext;
}


void _MeshNameHandler::printNames()
{
  printf("Input file: \n");
  printf("meshfile:    %s\n",_in_meshfile.c_str());
  printf("basename:    %s\n",_in_basename.c_str());
  printf("suffix:      %s\n",_in_extension.c_str());
  printf("path:        %s\n",_in_path.c_str());
  printf("Output file:\n");
}
