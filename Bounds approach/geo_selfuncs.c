/*-------------------------------------------------------------------------
 *
 * geo_selfuncs.c
 *	  Selectivity routines registered in the operator catalog in the
 *	  "oprrest" and "oprjoin" attributes.
 *
 * Portions Copyright (c) 1996-2020, PostgreSQL Global Development Group
 * Portions Copyright (c) 1994, Regents of the University of California
 *
 *
 * IDENTIFICATION
 *	  src/backend/utils/adt/geo_selfuncs.c
 *
 *	XXX These are totally bogus.  Perhaps someone will make them do
 *	something reasonable, someday.
 *
 *-------------------------------------------------------------------------
 */
#include "postgres.h"

#include <math.h>

#include "utils/builtins.h"
#include "utils/geo_decls.h"
#include "access/htup_details.h"
#include "catalog/pg_statistic.h"
#include "nodes/pg_list.h"
#include "optimizer/pathnode.h"
#include "optimizer/optimizer.h"
#include "utils/lsyscache.h"
#include "utils/typcache.h"
#include "utils/selfuncs.h"
#include "utils/rangetypes.h"

#include "catalog/pg_operator.h"
#include "catalog/pg_type.h"
#include "utils/float.h"
#include "utils/fmgrprotos.h"

static double calc_hist_selectivity_scalar(TypeCacheEntry *typcache,
										   const RangeBound *constbound,
										   const RangeBound *hist, int hist_nvalues,
										   bool equal);
static int	rbound_bsearch(TypeCacheEntry *typcache, const RangeBound *value,
						   const RangeBound *hist, int hist_length, bool equal);
static float8 get_position(TypeCacheEntry *typcache, const RangeBound *value,
						   const RangeBound *hist1, const RangeBound *hist2);
static float8 get_len_position(double value, double hist1, double hist2);
static float8 get_distance(TypeCacheEntry *typcache, const RangeBound *bound1,
						   const RangeBound *bound2);
static int	length_hist_bsearch(Datum *length_hist_values,
								int length_hist_nvalues, double value, bool equal);
static double calc_length_hist_frac(Datum *length_hist_values,
									int length_hist_nvalues, double length1, double length2, bool equal);

/*
 *	Selectivity functions for geometric operators.  These are bogus -- unless
 *	we know the actual key distribution in the index, we can't make a good
 *	prediction of the selectivity of these operators.
 *
 *	Note: the values used here may look unreasonably small.  Perhaps they
 *	are.  For now, we want to make sure that the optimizer will make use
 *	of a geometric index if one is available, so the selectivity had better
 *	be fairly small.
 *
 *	In general, GiST needs to search multiple subtrees in order to guarantee
 *	that all occurrences of the same key have been found.  Because of this,
 *	the estimated cost for scanning the index ought to be higher than the
 *	output selectivity would indicate.  gistcostestimate(), over in selfuncs.c,
 *	ought to be adjusted accordingly --- but until we can generate somewhat
 *	realistic numbers here, it hardly matters...
 */


/*
 * Selectivity for operators that depend on area, such as "overlap".
 */

Datum
areasel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.005);
}

Datum    //TODO
areajoinsel(PG_FUNCTION_ARGS)
{
    PlannerInfo *root = (PlannerInfo *) PG_GETARG_POINTER(0);
    Oid         operator = PG_GETARG_OID(1);
    List       *args = (List *) PG_GETARG_POINTER(2);
    JoinType    jointype = (JoinType) PG_GETARG_INT16(3);
    SpecialJoinInfo *sjinfo = (SpecialJoinInfo *) PG_GETARG_POINTER(4);
    Oid         collation = PG_GET_COLLATION();

    double      selec = 0.005;

    VariableStatData vardata1;
    VariableStatData vardata2;
    Oid         opfuncoid;
    AttStatsSlot sslot1;
	AttStatsSlot sslot2;
    int         nhist1;
	int         nhist2;
    RangeBound *hist_lower1;
    RangeBound *hist_upper1;
	RangeBound *hist_lower2;
    RangeBound *hist_upper2;
    int         i;
    Form_pg_statistic stats1 = NULL;
    Form_pg_statistic stats2 = NULL;
    TypeCacheEntry *typcache1 = NULL;
	TypeCacheEntry *typcache2 = NULL;
    bool        join_is_reversed;
    bool        empty;
	/*
	double quantity1;
	double quantity2;
	*/
	RangeType  *constrange = NULL;

    get_join_variables(root, args, sjinfo,
                       &vardata1, &vardata2, &join_is_reversed);

    typcache1 = range_get_typcache(fcinfo, vardata1.vartype);
	typcache2 = range_get_typcache(fcinfo, vardata2.vartype);
    opfuncoid = get_opcode(operator);

    memset(&sslot1, 0, sizeof(sslot1));
	memset(&sslot2, 0, sizeof(sslot2));

    /* Can't use the histogram with insecure range support functions */
    if (!statistic_proc_security_check(&vardata1, opfuncoid))
        PG_RETURN_FLOAT8((float8) selec);
    if (!statistic_proc_security_check(&vardata2, opfuncoid))
        PG_RETURN_FLOAT8((float8) selec);

    if (HeapTupleIsValid(vardata1.statsTuple))
    {
        stats1 = (Form_pg_statistic) GETSTRUCT(vardata1.statsTuple);/*
		quantity1 = stats1->stadistinct;
		printf("quantity2-1:%lf\n",quantity1);*/
        /* Try to get fraction of empty ranges */
        if (!get_attstatsslot(&sslot1, vardata1.statsTuple,
                             STATISTIC_KIND_BOUNDS_HISTOGRAM,
                             InvalidOid, ATTSTATSSLOT_VALUES))
        {
            ReleaseVariableStats(vardata1);
            ReleaseVariableStats(vardata2);
            PG_RETURN_FLOAT8((float8) selec);
        }
    }
	if (HeapTupleIsValid(vardata2.statsTuple))
    {
        stats2 = (Form_pg_statistic) GETSTRUCT(vardata2.statsTuple);/*
		quantity2 = stats2->stadistinct;
		printf("quantity2-2:%lf\n",quantity2);
		fflush(stdout);*/
        /* Try to get fraction of empty ranges */
        if (!get_attstatsslot(&sslot2, vardata2.statsTuple,
                             STATISTIC_KIND_BOUNDS_HISTOGRAM,
                             InvalidOid, ATTSTATSSLOT_VALUES))
        {
            ReleaseVariableStats(vardata1);
            ReleaseVariableStats(vardata2);
            PG_RETURN_FLOAT8((float8) selec);
        }
    }

    nhist1 = sslot1.nvalues;
	nhist2 = sslot2.nvalues;
    hist_lower1 = (RangeBound *) palloc(sizeof(RangeBound) * nhist1);
    hist_upper1 = (RangeBound *) palloc(sizeof(RangeBound) * nhist1);
	hist_lower2 = (RangeBound *) palloc(sizeof(RangeBound) * nhist2);
    hist_upper2 = (RangeBound *) palloc(sizeof(RangeBound) * nhist2);
    for (i = 0; i < nhist1; i++)
    {
        range_deserialize(typcache1, DatumGetRangeTypeP(sslot1.values[i]),
                          &hist_lower1[i], &hist_upper1[i], &empty);
        /* The histogram should not contain any empty ranges */
        if (empty)
            elog(ERROR, "bounds histogram contains an empty range");
    }
    for (i = 0; i < nhist2; i++)
    {
        range_deserialize(typcache2, DatumGetRangeTypeP(sslot2.values[i]),
                          &hist_lower2[i], &hist_upper2[i], &empty);
        /* The histogram should not contain any empty ranges */
        if (empty)
            elog(ERROR, "bounds histogram contains an empty range");
    }
	
	double bin_frac = 1.0/nhist1;	
	selec = 1.0;
	printf("bin_frac:%lf\n",bin_frac);
	
	for (i = 0; i < nhist1; i++)
    {
		double a,b;
		
		RangeBound bin_lower,bin_upper;
		bin_lower.inclusive = true;
		bin_lower.val = (Datum)(hist_lower1[i].val);
		bin_lower.infinite = false;
		bin_lower.lower = true;
		bin_upper.inclusive = true;
		bin_upper.val = (Datum)(hist_upper1[i].val);
		bin_upper.infinite = false;
		bin_upper.lower = false;

		//printf("%d a:%d  b:%d \n",i,bin_lower.val,bin_upper.val);

		a = calc_hist_selectivity_scalar(typcache2, &bin_lower, hist_upper2, nhist2, false);
		b = 1 - calc_hist_selectivity_scalar(typcache2, &bin_upper, hist_lower2, nhist2, true); 
		
		//printf("a:%lf  b:%lf \n",a,b);

		selec -= (a+b)*bin_frac;
		
    }
	printf("selec: %lf",selec);
		/*
    printf("hist_lower1 = [");
    for (i = 0; i < nhist1; i++)
    {
        printf("%d", DatumGetInt16(hist_lower1[i].val));
        if (i < nhist1 - 1)
            printf(", ");
    }
    printf("]\n");
    printf("hist_upper1 = [");
    for (i = 0; i < nhist1; i++)
    {
        printf("%d", DatumGetInt16(hist_upper1[i].val));
        if (i < nhist1 - 1)
            printf(", ");
    }
    printf("]\n");
	
	printf("hist_lower2 = [");
    for (i = 0; i < nhist2; i++)
    {
        printf("%d", DatumGetInt16(hist_lower2[i].val));
        if (i < nhist2 - 1)
            printf(", ");
    }
    printf("]\n");
    printf("hist_upper2 = [");
    for (i = 0; i < nhist2; i++)
    {
        printf("%d", DatumGetInt16(hist_upper2[i].val));
        if (i < nhist2 - 1)
            printf(", ");
    }
    printf("]\n");*/

    fflush(stdout);

    pfree(hist_lower1);
    pfree(hist_upper1);
    pfree(hist_lower2);
    pfree(hist_upper2);

    free_attstatsslot(&sslot1);    
	free_attstatsslot(&sslot2);

    ReleaseVariableStats(vardata1);
    ReleaseVariableStats(vardata2);

    CLAMP_PROBABILITY(selec);
    PG_RETURN_FLOAT8((float8) selec);
}

/*
 *	positionsel
 *
 * How likely is a box to be strictly left of (right of, above, below)
 * a given box?
 */

Datum
positionsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.1);
}

Datum
positionjoinsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.1);
}

/*
 *	contsel -- How likely is a box to contain (be contained by) a given box?
 *
 * This is a tighter constraint than "overlap", so produce a smaller
 * estimate than areasel does.
 */

Datum
contsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.001);
}

Datum
contjoinsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.001);
}



/*
 * Look up the fraction of values less than (or equal, if 'equal' argument
 * is true) a given const in a histogram of range bounds.
 */
static double
calc_hist_selectivity_scalar(TypeCacheEntry *typcache, const RangeBound *constbound,
							 const RangeBound *hist, int hist_nvalues, bool equal)
{
	Selectivity selec;
	int			index;

	/*
	 * Find the histogram bin the given constant falls into. Estimate
	 * selectivity as the number of preceding whole bins.
	 */
	index = rbound_bsearch(typcache, constbound, hist, hist_nvalues, equal);
	selec = (Selectivity) (Max(index, 0)) / (Selectivity) (hist_nvalues - 1);

	/* Adjust using linear interpolation within the bin */
	if (index >= 0 && index < hist_nvalues - 1)
		selec += get_position(typcache, constbound, &hist[index],
							  &hist[index + 1]) / (Selectivity) (hist_nvalues - 1);

	return selec;
}

/*
 * Binary search on an array of range bounds. Returns greatest index of range
 * bound in array which is less(less or equal) than given range bound. If all
 * range bounds in array are greater or equal(greater) than given range bound,
 * return -1. When "equal" flag is set conditions in brackets are used.
 *
 * This function is used in scalar operator selectivity estimation. Another
 * goal of this function is to find a histogram bin where to stop
 * interpolation of portion of bounds which are less than or equal to given bound.
 */
static int
rbound_bsearch(TypeCacheEntry *typcache, const RangeBound *value, const RangeBound *hist,
			   int hist_length, bool equal)
{
	int			lower = -1,
				upper = hist_length - 1,
				cmp,
				middle;

	while (lower < upper)
	{
		middle = (lower + upper + 1) / 2;
		cmp = range_cmp_bounds(typcache, &hist[middle], value);

		if (cmp < 0 || (equal && cmp == 0))
			lower = middle;
		else
			upper = middle - 1;
	}
	return lower;
}


/*
 * Binary search on length histogram. Returns greatest index of range length in
 * histogram which is less than (less than or equal) the given length value. If
 * all lengths in the histogram are greater than (greater than or equal) the
 * given length, returns -1.
 */
static int
length_hist_bsearch(Datum *length_hist_values, int length_hist_nvalues,
					double value, bool equal)
{
	int			lower = -1,
				upper = length_hist_nvalues - 1,
				middle;

	while (lower < upper)
	{
		double		middleval;

		middle = (lower + upper + 1) / 2;

		middleval = DatumGetFloat8(length_hist_values[middle]);
		if (middleval < value || (equal && middleval <= value))
			lower = middle;
		else
			upper = middle - 1;
	}
	return lower;
}

/*
 * Get relative position of value in histogram bin in [0,1] range.
 */
static float8
get_position(TypeCacheEntry *typcache, const RangeBound *value, const RangeBound *hist1,
			 const RangeBound *hist2)
{
	bool		has_subdiff = OidIsValid(typcache->rng_subdiff_finfo.fn_oid);
	float8		position;

	if (!hist1->infinite && !hist2->infinite)
	{
		float8		bin_width;

		/*
		 * Both bounds are finite. Assuming the subtype's comparison function
		 * works sanely, the value must be finite, too, because it lies
		 * somewhere between the bounds.  If it doesn't, arbitrarily return
		 * 0.5.
		 */
		if (value->infinite)
			return 0.5;

		/* Can't interpolate without subdiff function */
		if (!has_subdiff)
			return 0.5;

		/* Calculate relative position using subdiff function. */
		bin_width = DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
													 typcache->rng_collation,
													 hist2->val,
													 hist1->val));
		if (isnan(bin_width) || bin_width <= 0.0)
			return 0.5;			/* punt for NaN or zero-width bin */

		position = DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
													typcache->rng_collation,
													value->val,
													hist1->val))
			/ bin_width;

		if (isnan(position))
			return 0.5;			/* punt for NaN from subdiff, Inf/Inf, etc */

		/* Relative position must be in [0,1] range */
		position = Max(position, 0.0);
		position = Min(position, 1.0);
		return position;
	}
	else if (hist1->infinite && !hist2->infinite)
	{
		/*
		 * Lower bin boundary is -infinite, upper is finite. If the value is
		 * -infinite, return 0.0 to indicate it's equal to the lower bound.
		 * Otherwise return 1.0 to indicate it's infinitely far from the lower
		 * bound.
		 */
		return ((value->infinite && value->lower) ? 0.0 : 1.0);
	}
	else if (!hist1->infinite && hist2->infinite)
	{
		/* same as above, but in reverse */
		return ((value->infinite && !value->lower) ? 1.0 : 0.0);
	}
	else
	{
		/*
		 * If both bin boundaries are infinite, they should be equal to each
		 * other, and the value should also be infinite and equal to both
		 * bounds. (But don't Assert that, to avoid crashing if a user creates
		 * a datatype with a broken comparison function).
		 *
		 * Assume the value to lie in the middle of the infinite bounds.
		 */
		return 0.5;
	}
}


/*
 * Get relative position of value in a length histogram bin in [0,1] range.
 */
static double
get_len_position(double value, double hist1, double hist2)
{
	if (!isinf(hist1) && !isinf(hist2))
	{
		/*
		 * Both bounds are finite. The value should be finite too, because it
		 * lies somewhere between the bounds. If it doesn't, just return
		 * something.
		 */
		if (isinf(value))
			return 0.5;

		return 1.0 - (hist2 - value) / (hist2 - hist1);
	}
	else if (isinf(hist1) && !isinf(hist2))
	{
		/*
		 * Lower bin boundary is -infinite, upper is finite. Return 1.0 to
		 * indicate the value is infinitely far from the lower bound.
		 */
		return 1.0;
	}
	else if (isinf(hist1) && isinf(hist2))
	{
		/* same as above, but in reverse */
		return 0.0;
	}
	else
	{
		/*
		 * If both bin boundaries are infinite, they should be equal to each
		 * other, and the value should also be infinite and equal to both
		 * bounds. (But don't Assert that, to avoid crashing unnecessarily if
		 * the caller messes up)
		 *
		 * Assume the value to lie in the middle of the infinite bounds.
		 */
		return 0.5;
	}
}

/*
 * Measure distance between two range bounds.
 */
static float8
get_distance(TypeCacheEntry *typcache, const RangeBound *bound1, const RangeBound *bound2)
{
	bool		has_subdiff = OidIsValid(typcache->rng_subdiff_finfo.fn_oid);

	if (!bound1->infinite && !bound2->infinite)
	{
		/*
		 * Neither bound is infinite, use subdiff function or return default
		 * value of 1.0 if no subdiff is available.
		 */
		if (has_subdiff)
		{
			float8		res;

			res = DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
												   typcache->rng_collation,
												   bound2->val,
												   bound1->val));
			/* Reject possible NaN result, also negative result */
			if (isnan(res) || res < 0.0)
				return 1.0;
			else
				return res;
		}
		else
			return 1.0;
	}
	else if (bound1->infinite && bound2->infinite)
	{
		/* Both bounds are infinite */
		if (bound1->lower == bound2->lower)
			return 0.0;
		else
			return get_float8_infinity();
	}
	else
	{
		/* One bound is infinite, the other is not */
		return get_float8_infinity();
	}
}

/*
 * Calculate the average of function P(x), in the interval [length1, length2],
 * where P(x) is the fraction of tuples with length < x (or length <= x if
 * 'equal' is true).
 */
static double
calc_length_hist_frac(Datum *length_hist_values, int length_hist_nvalues,
					  double length1, double length2, bool equal)
{
	double		frac;
	double		A,
				B,
				PA,
				PB;
	double		pos;
	int			i;
	double		area;

	Assert(length2 >= length1);

	if (length2 < 0.0)
		return 0.0;				/* shouldn't happen, but doesn't hurt to check */

	/* All lengths in the table are <= infinite. */
	if (isinf(length2) && equal)
		return 1.0;

	/*----------
	 * The average of a function between A and B can be calculated by the
	 * formula:
	 *
	 *			B
	 *	  1		/
	 * -------	| P(x)dx
	 *	B - A	/
	 *			A
	 *
	 * The geometrical interpretation of the integral is the area under the
	 * graph of P(x). P(x) is defined by the length histogram. We calculate
	 * the area in a piecewise fashion, iterating through the length histogram
	 * bins. Each bin is a trapezoid:
	 *
	 *		 P(x2)
	 *		  /|
	 *		 / |
	 * P(x1)/  |
	 *	   |   |
	 *	   |   |
	 *	---+---+--
	 *	   x1  x2
	 *
	 * where x1 and x2 are the boundaries of the current histogram, and P(x1)
	 * and P(x1) are the cumulative fraction of tuples at the boundaries.
	 *
	 * The area of each trapezoid is 1/2 * (P(x2) + P(x1)) * (x2 - x1)
	 *
	 * The first bin contains the lower bound passed by the caller, so we
	 * use linear interpolation between the previous and next histogram bin
	 * boundary to calculate P(x1). Likewise for the last bin: we use linear
	 * interpolation to calculate P(x2). For the bins in between, x1 and x2
	 * lie on histogram bin boundaries, so P(x1) and P(x2) are simply:
	 * P(x1) =	  (bin index) / (number of bins)
	 * P(x2) = (bin index + 1 / (number of bins)
	 */

	/* First bin, the one that contains lower bound */
	i = length_hist_bsearch(length_hist_values, length_hist_nvalues, length1, equal);
	if (i >= length_hist_nvalues - 1)
		return 1.0;

	if (i < 0)
	{
		i = 0;
		pos = 0.0;
	}
	else
	{
		/* interpolate length1's position in the bin */
		pos = get_len_position(length1,
							   DatumGetFloat8(length_hist_values[i]),
							   DatumGetFloat8(length_hist_values[i + 1]));
	}
	PB = (((double) i) + pos) / (double) (length_hist_nvalues - 1);
	B = length1;

	/*
	 * In the degenerate case that length1 == length2, simply return
	 * P(length1). This is not merely an optimization: if length1 == length2,
	 * we'd divide by zero later on.
	 */
	if (length2 == length1)
		return PB;

	/*
	 * Loop through all the bins, until we hit the last bin, the one that
	 * contains the upper bound. (if lower and upper bounds are in the same
	 * bin, this falls out immediately)
	 */
	area = 0.0;
	for (; i < length_hist_nvalues - 1; i++)
	{
		double		bin_upper = DatumGetFloat8(length_hist_values[i + 1]);

		/* check if we've reached the last bin */
		if (!(bin_upper < length2 || (equal && bin_upper <= length2)))
			break;

		/* the upper bound of previous bin is the lower bound of this bin */
		A = B;
		PA = PB;

		B = bin_upper;
		PB = (double) i / (double) (length_hist_nvalues - 1);

		/*
		 * Add the area of this trapezoid to the total. The point of the
		 * if-check is to avoid NaN, in the corner case that PA == PB == 0,
		 * and B - A == Inf. The area of a zero-height trapezoid (PA == PB ==
		 * 0) is zero, regardless of the width (B - A).
		 */
		if (PA > 0 || PB > 0)
			area += 0.5 * (PB + PA) * (B - A);
	}

	/* Last bin */
	A = B;
	PA = PB;

	B = length2;				/* last bin ends at the query upper bound */
	if (i >= length_hist_nvalues - 1)
		pos = 0.0;
	else
	{
		if (DatumGetFloat8(length_hist_values[i]) == DatumGetFloat8(length_hist_values[i + 1]))
			pos = 0.0;
		else
			pos = get_len_position(length2, DatumGetFloat8(length_hist_values[i]), DatumGetFloat8(length_hist_values[i + 1]));
	}
	PB = (((double) i) + pos) / (double) (length_hist_nvalues - 1);

	if (PA > 0 || PB > 0)
		area += 0.5 * (PB + PA) * (B - A);

	/*
	 * Ok, we have calculated the area, ie. the integral. Divide by width to
	 * get the requested average.
	 *
	 * Avoid NaN arising from infinite / infinite. This happens at least if
	 * length2 is infinite. It's not clear what the correct value would be in
	 * that case, so 0.5 seems as good as any value.
	 */
	if (isinf(area) && isinf(length2))
		frac = 0.5;
	else
		frac = area / (length2 - length1);

	return frac;
}

