#ifndef FILE_OPERATION_H
///Define this macro to prevent from including this header file more than once.
#define FILE_OPERATION_H

#include "stl.h"
#include "type.h"

inline bool is_directory(const string directory) {
	if (directory.length() == 0) return false;
#ifdef WINDOWS
	return directory[directory.length() - 1] == '\\';
#else
	return directory[directory.length() - 1] == '/';
#endif
}

inline bool file_exists (const string file_name) {
  ifstream fin;
  // We must use the string class's c_str() function to convert
  // The C++ string into an old C string (character array).
  // This is because the open() function expects a character array
  fin.open (file_name.c_str());
  if (fin.fail()) return false;
  fin.close();
  return true;
}

inline bool complete_directory(string &directory) {
	int l = (int)directory.length();
#ifdef WINDOWS
	if (directory.find("/") != string::npos) return false;
	if (directory[l-1] == '\\') return true;
	else {
		directory = directory + "\\";
		return true;
	}
#else
	if (directory.find("\\") != string::npos) return false;
	if (directory[l-1] == '/') return true;
	else {
		directory = directory + "/";
		return true;
	}
#endif
}

inline string get_path(const string full_pathname){
	size_t i = full_pathname.find_last_of("\\");
	if (i == string::npos || i == full_pathname.length() - 1) return full_pathname; else return full_pathname.substr(0, i+1);
}

inline string get_file_name(const string file_name) {
	size_t i = file_name.find_last_of("\\");
	if (i == string::npos) return file_name;
	else if (i == file_name.length() - 1) return "";
	else return file_name.substr(i+1);
}

inline bool load_file(const string file_name, char* &buf, long &filesize){
	FILE *in;

	if (!file_exists(file_name)) return false;
	in = fopen(file_name.c_str(),"rb");  // read binary
	if ( !in ) return  false;
		// determine file size
	fseek(in,0,SEEK_END);
	filesize = ftell(in);
	fseek(in,0,SEEK_SET);

	// allocate buffer and read in file contents
	buf = new char[filesize+1];
	size_t rl = fread(buf,sizeof(char),filesize,in);
	buf[filesize] = 0;
	fclose(in);
	return true;
}

inline bool write_file(const string file_name, const char* const buf, long size) {
	FILE *out;

	out = fopen(file_name.c_str(),"wb");  // read binary
	if ( !out ) return  false;

	fwrite(buf,sizeof(char),size,out);
	fclose(out);
	return true;
}

#ifdef WINDOWS

//#include "shared_func.h"
#include <windows.h>

inline BOOL get_file_time(const string filename, LPFILETIME lpCreationTime, LPFILETIME lpLastAccessTime, LPFILETIME lpLastWriteTime){
	HANDLE hFile; 
    hFile = CreateFile(filename.c_str(),               // file to open
                       GENERIC_READ,          // open for reading
                       FILE_SHARE_READ,       // share for reading
                       NULL,                  // default security
                       OPEN_EXISTING,         // existing file only
                       FILE_ATTRIBUTE_NORMAL, // normal file
                       NULL);                 // no attr. template
 
    if (hFile == INVALID_HANDLE_VALUE) 
    { 
        printf("Could not open file (error %d)\n", GetLastError());
        return false; 
    }
	BOOL result = GetFileTime(hFile, lpCreationTime, lpLastAccessTime, lpLastWriteTime);
	if (!result) {
        printf("Could not get file time (error %d)\n", GetLastError());
	}
    CloseHandle(hFile);
	return result;
}

inline FILETIME get_file_modify_time(const string file) {
	FILETIME filetime1, filetime2, filetime3;
	get_file_time(file, &filetime1, &filetime2, &filetime3);
	return filetime3;
}

inline bool file_changed(const string file, const FILETIME filetime) {
	FILETIME filetime1, filetime2, filetime3;
	get_file_time(file, &filetime1, &filetime2, &filetime3);
	return (0!=CompareFileTime(&filetime, &filetime3));
}

inline LONG compare_file_time(const string file1, const string file2) {
	FILETIME filetime11, filetime12, filetime13, filetime21, filetime22, filetime23;
	get_file_time(file1, &filetime11, &filetime12, &filetime13);
	get_file_time(file2, &filetime21, &filetime22, &filetime23);
	return CompareFileTime(&filetime13, &filetime23);
}

inline BOOL get_file_size(const string filename, int64 *file_size){
	HANDLE hFile; 
    hFile = CreateFile(filename.c_str(),               // file to open
                       GENERIC_READ,          // open for reading
                       FILE_SHARE_READ,       // share for reading
                       NULL,                  // default security
                       OPEN_EXISTING,         // existing file only
                       FILE_ATTRIBUTE_NORMAL, // normal file
                       NULL);                 // no attr. template
 
    if (hFile == INVALID_HANDLE_VALUE) 
    { 
        printf("Could not open file (error %d)\n", GetLastError());
        return false; 
    }
	BOOL result = GetFileSizeEx(hFile, (PLARGE_INTEGER)file_size);
	if (!result) {
        printf("Could not get file size (error %d)\n", GetLastError());
	}
    CloseHandle(hFile);
	return result;
}

#endif //#ifdef WINDOWS

#endif // FILE_OPERATION_H
