
 typedef struct{
  double prob[5];
/*--------------------probabilities:
  array-style indexing: if (i,j) = current node then 
  stays:   prob[0]      i,j same   
  west:    prob[1]      j decreases by 1
  east:    prob[2]      j increases by 1
  north:   prob[3]      i decreases by 1
  south:   prob[4]      i increases by 1
*/
} PointProb;

typedef struct _dom{
  int ni;       
  int nj;
  int qi;
  int qj; 
  int westNB;
  int eastNB;
  int southNB;
  int northNB;
  double* west;
  double* east;
  double* south;
  double* north;
} myDomain, *DomainPtr;
