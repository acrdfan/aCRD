// Diffusivity Library
// In MOD_LIB_aCRD
// Version 0 alpha demo
// J. Fan 2018

#ifdef DIFFUSIVITYLIB_EXPORTS
#define DIFFUSIVITYLIB_API __declspec(dllexport)
#else
#define DIFFUSIVITYLIB_API __declspec(dllimport)
#endif

#ifdef __cplusplus
extern "C" {
#endif
	DIFFUSIVITYLIB_API double getDvValue(double E, int type, long size);
	extern DIFFUSIVITYLIB_API double nDvLibVersion;
#ifdef __cplusplus
}
#endif
