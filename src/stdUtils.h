#ifndef StreamUtils_h__
#define StreamUtils_h__
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <istream>
#include <iostream>
#include <iomanip>
#include <algorithm>

static const double DEG2RAD = 0.0174532925199433;
static const double RAD2DEG = 57.2957795130823;
static const double PI = 3.14159265358979;


#ifdef WIN32
    #include <windows.h>
    #include <shlobj.h>
#else
    #include  <sys/stat.h>
    #include <dirent.h>
#endif
// trim from left
static inline std::string &ltrim(std::string &s, const char* t = " \t\n\r\f\v")
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}
// trim from right
static inline std::string &rtrim(std::string &s, const char* t = " \t\n\r\f\v")
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}
// trim from left & right
static inline std::string &trim(std::string &s, const char* t = " \t\n\r\f\v")
{
    return ltrim(rtrim(s, t), t);
}

//http://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
template <typename T>
static inline std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}
template <typename T>
static inline std::string STR(const T value)
{
	return std::to_string(value);
}
template <typename T>
static inline std::string STR(const T value, const int n)
{
	return to_string_with_precision(value);
}
template <typename T>
static inline std::string STR31(const T * pvalue, const int precision = 6, const char* delimiter = " ")
{
	std::ostringstream out;
	out.precision(precision);
	out << pvalue[0] << delimiter << pvalue[1] << delimiter << pvalue[2];
	return out.str();
}
static inline std::vector<std::vector<std::string> > ReadSpaceDelimitedFile(const char* file)
{
    std::vector<std::vector<std::string> > data;
    std::ifstream infile(file);
    if (!infile.is_open())
    {
        std::cout << "could not open the file " << file << std::endl;
        return data;
    }
    while (infile)
    {
        std::string s;
        if (!getline(infile, s)) break;
        
        std::istringstream ss(s);
        std::vector <std::string> record;
        
        while (ss)
        {
            std::string s;
            
            ss >> s;
            
            //if (!getline(ss, s, ',')) break;
            if (!trim(s).empty())
                record.push_back(s);
        }
        
        data.push_back(record);
    }
    if (!infile.eof())
    {
        std::cerr << "Fooey!\n";
    }
    
    return data;
}
static inline std::string GetFileName(const std::string &FilePath)
{
    size_t found;
    found = FilePath.find_last_of("/\\");
    return FilePath.substr(found + 1);
}

static inline std::string GetFileNameWithoutExtention(const std::string &FilePath)
{
    std::string filename = GetFileName(FilePath);
    size_t found;
    found = filename.find_last_of(".");
    return filename.substr(0,found);
}
static inline std::string GetFileNameExtention(const std::string &FilePath)
{
    std::string filename = GetFileName(FilePath);
    size_t found;
    found = filename.find_last_of(".");
    return filename.substr(found+1,filename.size());
}
static inline std::string GetFileFolderName(const std::string &FilePath)
{
    size_t found;
    found = FilePath.find_last_of("/\\");
    
    return FilePath.substr(0,found);
}
//returns the full tile path, but just drop the extension
static inline std::string GetFilePathWithoutExtention(const std::string &FilePath)
{

  size_t found;
  found = FilePath.find_last_of(".");
  return FilePath.substr(0, found);
}


//https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
static inline std::vector<std::string> split(const std::string& s, const char delimiter=' ')
{
	std::vector<std::string> tokens;

	for (size_t start = 0, end; start < s.length(); start = end + 1)
	{
		size_t position = s.find(delimiter, start);
		end = position != std::string::npos ? position : s.length();
        std::string sb=s.substr(start, end - start);
        std::string& token =trim(sb);
		if (!token.empty())
		{
			tokens.push_back(token);
		}
	}
	return tokens;
}
//https://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int
static inline bool has_only_digits(const std::string &s) {return s.find_first_not_of("0123456789") == std::string::npos;}

static inline std::string GetMyDocsDirA()
{
    std::string my_docs_dir;
    
    #ifdef WIN32
    char szPath[MAX_PATH];
    if (SUCCEEDED(SHGetFolderPathA(NULL,
                                   CSIDL_PERSONAL | CSIDL_FLAG_CREATE,
                                   NULL,
                                   0,
                                   szPath)))
    {
        my_docs_dir = szPath;
    }
    #endif
    return my_docs_dir;
};
static inline std::string GetExePathA()
{
  #ifdef WIN32
  char szPath[_MAX_PATH + 1];
  GetModuleFileNameA(0, szPath, _MAX_PATH + 1);
  
  return std::string(szPath);
  #else
  char arg1[20];
  char exepath[PATH_MAX + 1] = { 0 };
  
  sprintf(arg1, "/proc/%d/exe", getpid());
  readlink(arg1, exepath, 1024);
  return std::string(exepath);;
  #endif
}
inline std::string ReadWholeFile(const std::string& FilePath)
{
    //read the whole file
    std::ifstream ifs(FilePath.c_str());
    ifs.seekg(0, std::ios::end);
    size_t length = ifs.tellg();
    std::string Buffer(length, '0');
    ifs.seekg(0, std::ios::beg);
    ifs.read(&Buffer[0], length);
    ifs.close();
    
    return Buffer;
    
}
static inline std::vector<std::string> GetFilesInDirectoryA(const std::string directory, const std::string filter="*")
{
    std::vector<std::string> out;
    #ifdef WIN32
    HANDLE dir;
    WIN32_FIND_DATAA file_data;
    
    if ((dir = FindFirstFileA((directory + "/" + filter).c_str(), &file_data)) == INVALID_HANDLE_VALUE)
        return out; /* No files found */
        
    do
    {
        const std::string file_name = file_data.cFileName;
        const std::string full_file_name = directory + "/" + file_name;
        const bool is_directory = (file_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) != 0;
        
        if (file_name[0] == '.')
            continue;
            
        if (is_directory)
            continue;
            
        out.push_back(full_file_name);
    }
    while (FindNextFileA(dir, &file_data));
    
    FindClose(dir);
    
    return out;
    #else
    DIR* dir;
    class dirent* ent;
    class stat st;
    
    dir = opendir(directory.c_str());
    while ((ent = readdir(dir)) != NULL)
    {
        const std::string file_name = ent->d_name;
        const std::string full_file_name = directory + "/" + file_name;
    
        if (file_name[0] == '.')
            continue;
    
        if (stat(full_file_name.c_str(), &st) == -1)
            continue;
    
        const bool is_directory = (st.st_mode & S_IFDIR) != 0;
    
        if (is_directory)
            continue;
    
        out.push_back(full_file_name);
    }
    closedir(dir);
    return out;
    #endif
} // GetFilesInDirectory


////////////////////////////////////////////////////////////////////////// VECTOR UTILS

template <typename T>
inline std::vector<size_t> sort_indexes_asc(const std::vector<T> &v)
{

    // initialize original index locations
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
    [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2]; });
    
    return idx;
}
template <typename T>
inline std::vector<size_t> sort_indexes_des(const std::vector<T> &v)
{

    // initialize original index locations
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
    [&v](std::size_t i1, std::size_t i2) {return v[i1] > v[i2]; });
    
    return idx;
}

//http://stackoverflow.com/questions/9323903/most-efficient-elegant-way-to-clip-a-number
template <typename T>
T inline clip(const T& n, const T& lower, const T& upper)
{
    return std::max(lower, std::min(n, upper));
}


//http://stackoverflow.com/questions/6892754/creating-a-simple-configuration-file-and-parser-in-c
inline std::map<std::string, std::string> ReadConfig(const std::string &ConfigFile)
{
    std::ifstream is_file(ConfigFile.c_str());
	std::map<std::string, std::string> config;
	if (!is_file.is_open())
    {
        std::cout << "could not find the config file " << ConfigFile << std::endl;
        return config;
    }   
    std::string line;
    while (std::getline(is_file, line))
    {
		if (line.find('#')!=std::string::npos)
        continue;
        std::istringstream is_line(line);
        std::string key;
        if (std::getline(is_line, key, '='))
        {
            key = trim(key);
            std::string value;
            if (std::getline(is_line, value))
            {
                if (config.count(key))
                    std::cout << "Warning duplicate setting " << key << std::endl;
                    
                config[key] = trim(value);
            }
            
        }
    }
    is_file.close();
    
    return config;
}

////////////////////////////////////////////////////////////////////////// MATH UTILS
template<typename T>
static inline T DotProduct31(const T* const A, const T* const B)
{
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
}
template<typename T>
static inline void Subtract31(T* AminusB,const T* const A, const T* const B)
{
    AminusB[0] = A[0] - B[0];
    AminusB[1] = A[1] - B[1];
    AminusB[2] = A[2] - B[2];
}
//A-=B
template<typename T>
static inline void Subtract31(T* AminusB, const T* const B)
{
	AminusB[0] -= B[0];
	AminusB[1] -= B[1];
	AminusB[2] -= B[2];
}
template<typename T, int n>
static inline void Subtractn(T* AminusB, const T* const A, const T* const B)
{
	for(int i=0;i<n;++i)
	AminusB[i] = A[i] - B[i];	
}

template<typename T>
static inline T Distance31(const T* const vec1, const T* const vec2)
{
    T res[3];
    Subtract31(res,vec1, vec2);
    return sqrt(DotProduct31(res, res));
}

template<typename T>
static inline void product441(T* result, const T* const mat44, const T* const mat41)
{
	result[0] = mat44[0]  * mat41[0] + mat44[1]  * mat41[1] + mat44[2]  * mat41[2] + mat44[3]  * mat41[3];
	result[1] = mat44[4]  * mat41[0] + mat44[5]  * mat41[1] + mat44[6]  * mat41[2] + mat44[7]  * mat41[3];
	result[2] = mat44[8]  * mat41[0] + mat44[9]  * mat41[1] + mat44[10] * mat41[2] + mat44[11] * mat41[3];
	result[3] = mat44[12] * mat41[0] + mat44[13] * mat41[1] + mat44[14] * mat41[2] + mat44[15] * mat41[3];
}
template<typename T>
static inline void product331(T* result, const T* const mat33, const T* const mat31)
{
	result[0] = mat33[0] * mat31[0] + mat33[1] * mat31[1] + mat33[2] * mat31[2]  ;
	result[1] = mat33[3] * mat31[0] + mat33[4] * mat31[1] + mat33[5] * mat31[2]  ;
	result[2] = mat33[6] * mat31[0] + mat33[7] * mat31[1] + mat33[8] * mat31[2] ;
}
template<typename T>
static inline void product33(T* result, const T* const left33, const T* const right33)
{
	result[0] = left33[0] * right33[0] + left33[1] * right33[3] + left33[2] * right33[6];
	result[1] = left33[0] * right33[1] + left33[1] * right33[4] + left33[2] * right33[7];
	result[2] = left33[0] * right33[2] + left33[1] * right33[5] + left33[2] * right33[8];

	result[3] = left33[3] * right33[0] + left33[4] * right33[3] + left33[5] * right33[6];
	result[4] = left33[3] * right33[1] + left33[4] * right33[4] + left33[5] * right33[7];
	result[5] = left33[3] * right33[2] + left33[4] * right33[5] + left33[5] * right33[8];

	result[6] = left33[6] * right33[0] + left33[7] * right33[3] + left33[8] * right33[6];
	result[7] = left33[6] * right33[1] + left33[7] * right33[4] + left33[8] * right33[7];
	result[8] = left33[6] * right33[2] + left33[7] * right33[5] + left33[8] * right33[8];
}


template<typename T>
static inline void ComputeMean(int n, const T* const Values, double MeanMinMax[3])
{
    double Sum = 0;
    T Max = -std::numeric_limits<T>::max();
    T Min = std::numeric_limits<T>::max();
    
    for (int i = 0; i < n; i++)
    {
        if (Values[i] > Max)
        {
            Max = Values[i];
        }
        else if (Values[i] < Min)
        {
            Min = Values[i];
        }
        
        Sum += Values[i];
    }
    
    MeanMinMax[0] = Sum / n;
    MeanMinMax[1] = Min;
    MeanMinMax[2] = Max;
    
}
template<typename T>
static inline double ComputeMean(int n, const T* const Values)
{
    double Sum = 0;
    for (int i = 0; i < n; ++i)
        Sum += Values[i];
        
    return Sum / n;
}
template<typename T>
static inline void Scale31(T scale, T * values)
{
	values[0] *= scale;
	values[1] *= scale;
	values[2] *= scale;
}
template<typename T>
static inline void Scale31(T* result, T scale, const  T * values)
{
	result[0]=values[0] * scale;
	result[1]=values[1] * scale;
	result[2]=values[2] * scale;
}

template<typename T>
static inline void Scale41(T scale, T * values)
{
	values[0] *= scale;
	values[1] *= scale;
	values[2] *= scale;
	values[3] *= scale;
}
template<typename T>
static inline void Scale9(T scale, T * values)
{
	values[0] *= scale;	values[1] *= scale;	values[2] *= scale;
	values[3] *= scale;	values[4] *= scale;	values[5] *= scale;
	values[6] *= scale;	values[7] *= scale;	values[8] *= scale;
}
template<typename T, int n>
static inline void Scale(T scale, T * values)
{	
	for (int i = 0; i < n; ++i)
		values[i] *= scale;
}

template<typename T>
static inline T Norm31(const T * const values)
{
	return sqrt(values[0] *values[0] + values[1] * values[1] + values[2] * values[2]);
}
template<typename T>
static inline T Norm41(const T * const values)
{
	return sqrt(values[0] * values[0] + values[1] * values[1] + values[2] * values[2]+ values[3] * values[3]);
}

//A=B
template<typename T>
static inline void Copy31(T * dst, const T* src)
{
	dst[0] = src[0];
	dst[1] = src[1];
	dst[2] = src[2];
}
template<typename T>
static inline void Copy41(T * dst, const T* src)
{
	dst[0] = src[0];
	dst[1] = src[1];
	dst[2] = src[2];
	dst[3] = src[3];
}

//A=A+alpha*B
template<typename T>
static inline void Sum31(T * A, const T* B, T alpha = T(1.0))
{
	A[0] += alpha*B[0];
	A[1] += alpha*B[1];
	A[2] += alpha*B[2];
}
//A=A+alpha*B
template<typename T>
static inline void Sum41(T * A, const T* B, T alpha = T(1.0))
{
	A[0] += alpha*B[0];
	A[1] += alpha*B[1];
	A[2] += alpha*B[2];
	A[3] += alpha*B[3];
}
//Result=A+alpha*B
template<typename T>
static inline void Sum31(T* Result,  const T * A, const T* B,  T alpha=T(1.0))
{
	Result[0] =A[0]+ alpha*B[0];
	Result[1] =A[1]+ alpha*B[1];
	Result[2] =A[2]+ alpha*B[2];

}
//Result=A+alpha*B
template<typename T>
static inline void Sum41(T* Result, const T * A, const T* B, T alpha = T(1.0))
{
	Result[0] = A[0] + alpha*B[0];
	Result[1] = A[1] + alpha*B[1];
	Result[2] = A[2] + alpha*B[2];
	Result[3] = A[3] + alpha*B[3];

}
template<typename T> static inline void setzero3(T * val) { val[0] = val[1] = val[2] = T(0.0); }
template<typename T> static inline void setzero4(T * val) { val[0] = val[1] = val[2] = val[3]= T(0.0); }

/*skew symmetric matrix*/
template<typename T>
static inline void skew(T * skewt, const T * t)
{
	skewt[0] = T(0.0);	skewt[1] = -t[2];	skewt[2] = t[1];
	skewt[3] = t[2]; 	skewt[4] = T(0.0); 		skewt[5] = -t[0];
	skewt[6] = -t[1]; 	skewt[7] = t[0]; 	    skewt[8] = T(0.0);
}
/*skew symmetric matrix squared*/
template<typename T>
static inline void skew2(T * sewsqt, const T * t)
{

	sewsqt[0] = -(t[2] * t[2] + t[1] * t[1]);		sewsqt[1] = t[1] * t[0];	sewsqt[2] = t[0] * t[2];
	sewsqt[3] = sewsqt[1]; 	sewsqt[4] = -(t[2] * t[2] + t[0] * t[0]); 		sewsqt[5] = t[2] * t[1];
	sewsqt[6] = sewsqt[2]; 	sewsqt[7] = sewsqt[5]; 	sewsqt[8] = -(t[1] * t[1] + t[0] * t[0]);
}

//http://stackoverflow.com/questions/8942950/how-do-i-find-the-orthogonal-projection-of-a-point-onto-a-plane
template<typename T>
static inline void PointToPlaneProjection(const T* const Plane4Coef, const T* const point, T Projection[3])
{
    // a vector from arbitrary point on the plane to the point
    T vec[3] = { 0 };
    
    // find the largest element of normal vector to the plane
    int largestElement = -1;
    if (abs(Plane4Coef[0]) >= abs(Plane4Coef[1]) && abs(Plane4Coef[0]) >= abs(Plane4Coef[2]))
    {
        largestElement = 0;
    }
    else if (abs(Plane4Coef[1]) >= abs(Plane4Coef[0]) && abs(Plane4Coef[1]) >= abs(Plane4Coef[2]))
    {
        largestElement = 1;
    }
    else if (abs(Plane4Coef[2]) >= abs(Plane4Coef[0]) && abs(Plane4Coef[2]) >= abs(Plane4Coef[1]))
    {
        largestElement = 2;
    }
    
    vec[0] = point[0]; vec[1] = point[1]; vec[2] = point[2];
    
    vec[largestElement] += Plane4Coef[3] / Plane4Coef[largestElement];
    
    // plane normal vector check to be unit
    T norm = sqrt(Plane4Coef[0] * Plane4Coef[0] + Plane4Coef[1] * Plane4Coef[1] + Plane4Coef[2] * Plane4Coef[2]);
    
    T dotProduct = (vec[0] * Plane4Coef[0] + vec[1] * Plane4Coef[1] + vec[2] * Plane4Coef[2]) / norm;
    
    
    Projection[0] = point[0] - dotProduct*Plane4Coef[0] / norm;
    Projection[1] = point[1] - dotProduct*Plane4Coef[1] / norm;
    Projection[2] = point[2] - dotProduct*Plane4Coef[2] / norm;
    
}
//https://math.stackexchange.com/questions/1300484/distance-between-line-and-a-point
template<typename T>
static inline void PointToLineProjection(const T *const LineDirection, const T *const PointOnline, const T *const pointOffline, T PointToLineVector[3])
{
  T dotvec = DotProduct31(LineDirection, LineDirection);
  
  
  Subtract31(PointToLineVector,pointOffline, PointOnline);
  T dot = DotProduct31(PointToLineVector, LineDirection);
  
  PointToLineVector[0] -=  dot / dotvec*PointToLineVector[0];
  PointToLineVector[1] -=  dot / dotvec*PointToLineVector[1];
  PointToLineVector[2] -=  dot / dotvec*PointToLineVector[2];
  
}

//https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
template<typename T>
static inline bool LineToPlaneIntersection(const T *const LineDirection, const T *const PointOnline, const T *const Plane4Coeef, T IntersectionPoint[3])
{

  T dotln = DotProduct31(LineDirection, Plane4Coeef);
  if (dotln == T(0.0))
    return false;
    
  T origin[3] = { 0.0 };
  T P0[3]; //a point on a plane, comes from projection of the origin into plane
  PointToPlaneProjection(Plane4Coeef, origin, P0);
  
  T vec[3];
  Subtract31(vec,P0, PointOnline);
  T dot= DotProduct31(vec, Plane4Coeef);
  
  T d = dot / dotln;
  
  IntersectionPoint[0] = d*LineDirection[0] + PointOnline[0];
  IntersectionPoint[1] = d*LineDirection[1] + PointOnline[1];
  IntersectionPoint[2] = d*LineDirection[2] + PointOnline[2];
  
  return true;
  
}
#endif // StreamUtils_h__
