#ifndef STRING_OPERATION_H
///Define this macro to prevent from including this header file more than once.
#define STRING_OPERATION_H

#include "stl.h"
#include "math_utils.h"

inline vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true) {
    // Skip delimiters at beginning.
	string::size_type lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);
	vector<string> result;
	result.clear();

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
		//__ASSERT(pos > lastPos || !skip_empty, "internal error, pos <= lastPos.\n");

		//if (pos == lastPos) result.push_back("");
		result.push_back(str.substr(lastPos, pos - lastPos));

		if (pos == string::npos) break;
		if (pos == str.length() - 1) {
			if (!skip_empty) result.push_back("");
			break;
		}
        // Skip delimiters.  Note the "not_of"
		lastPos = skip_empty ? str.find_first_not_of(delimiters, pos) : pos + 1;
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
	return result;
}

inline string tolower(const string s)
{
    string result = "";
	int i;
	for (i = 0; i < (int)s.length(); i++)
		result += (char)tolower(s[i]);
    return result;
}

inline string toupper(const string s)
{
    string result = "";
	int i;
	for (i = 0; i < (int)s.length(); i++)
		result += (char)toupper(s[i]);
    return result;
}

inline string substitute(const string s, const string from, const string to) {
	string p;
	if ( from.length() == 0 )
	{
		p = s;
	}
	else
	{
		int i;
		int istart;
	
		istart = 0;
		while (istart < (int)s.length())
		{
			i = (int)s.find(from, istart);
			if (i == (int)string::npos)
			{
				p += s.substr(istart);
				break;
			}
			p += s.substr(istart, i - istart);
			p += to;
			istart = i + (int)from.length();
		}
	}
	return p;
}

inline string int2str(const int i){
	char buf[1024];
	sprintf(buf, "%d", i);
	return string(buf);
}

inline int str2int(const string s){
	int i;
	sscanf(s.c_str(), "%d", &i);
	return i;
}

inline string int_vec2str(const vector<int> nums) {
	string result = "";
	for (int i = 0; i < (int)nums.size(); i++) {
		if (i > 0) result += ",";
		result += int2str(nums[i]);
	}
	return result;
}

inline vector<int> str2int_vec(const string s){
	vector<int> results;
	vector<string> nums = string_tokenize(s, ", ");
	for (int i = 0; i < (int)nums.size(); i++) {
		results.push_back(str2int(nums[i]));
	}
	return results;
}

inline bool is_int(const string s){
	return (s == int2str(str2int(s)));
}

inline string double2str(const double d){
	char buf[1024];
	sprintf(buf, "%lf", d);
	return string(buf);
}

inline double str2double(const string s){
	double d;
	sscanf(s.c_str(), "%lf", &d);
	return d;
}

inline bool is_num(const string s){
	return (strspn(s.c_str(), "01234567890eE.+-") == s.length() && strcspn(s.c_str(), "01234567890") < s.length());
}

inline string trim_space(const string s, const string& delimiters = " \t\n\r") {
	string::size_type pos1 = s.find_first_not_of(delimiters);
	string::size_type pos2 = s.find_last_not_of(delimiters);
	if (pos1 != string::npos && pos2 != string::npos && pos2 >= pos1) {
		return s.substr(pos1, pos2 - pos1 + 1);
	} else return "";
}

inline bool is_empty(const string s) {
	return (trim_space(s) == "");
}

inline string get_random_string(){
	return int2str(rand()*rand());
}

inline string get_time_string(){
	return int2str((int)time(NULL));
}

inline string current_time(){
	time_t rawtime;
	time(&rawtime);
	return string(ctime(&rawtime));
}

#endif // STRING_OPERATION_H
