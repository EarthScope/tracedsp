#ifndef IIRFILTER_H
#define IIRFILTER_H 1

#ifdef __cplusplus
extern "C" {
#endif

/* Recognized data types for input and output:
 * 'i' : 32-bit integers
 * 'f' : 32-bit floating point numbers
 * 'd' : 64-bit floating point numbers
 */

int
iirfilter ( void *input, char inputtype,
	    int samplecnt, int reverse,
	    void **output, char outputtype, 
	    int highorder, double highcutoff,
	    int loworder, double lowcutoff,
	    double samprate, int verbose );

#ifdef __cplusplus
}
#endif

#endif /* IIRFILTER_H */
