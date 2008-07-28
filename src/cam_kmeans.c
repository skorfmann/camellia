#include <stdlib.h>
#include <stdio.h>

#define NB_DIMENSIONS 4

typedef struct {
    int s[NB_DIMENSIONS];
} CamColorSample;

typedef struct {
    CamColorSample point;
    CamColorSample wgtCent;
    int count;
} CamColorCenter;

typedef struct _CamColorKdTreeNode {
    int dim;			// dimension 
    int value;			// edge value
    int count;			// number of samples in the cell
    CamColorSample wgtCent;	// weight of all samples in the cell
    struct _CamColorKdTreeNode *right; // direction to the right child in the kdtree (left child is + 1). NULL for leaves
    CamColorSample *point;
} CamColorKdTreeNode;

#define IS_LEAF(x) ((x)->right == NULL)

CamColorKdTreeNode *camColorKdTreeRecurse(CamColorKdTreeNode *parent, CamColorSample **samples, int nbs, int nb_dim)
{
    int i, dim, dsplit, leaf = 0, nbs_left = 0;
    CamColorKdTreeNode *retval;
    CamColorSample *p;
    int d, avg, bavg, avgdev, bavgdev;
    
    if (nbs > 1) { // If there is more than one sample concerned by this recursion
	// Compute the average and average deviation in all dimensions
	for (dim = 0; dim < nb_dim; dim++) {
	    avg = samples[0]->s[dim];
	    for (i = 1; i < nbs; i++) avg += samples[i]->s[dim];
	    avg /= nbs; // Warning : division
	    d = samples[0]->s[dim] - avg;
	    avgdev = d * d;
	    for (i = 1; i < nbs; i++) {
		d = samples[i]->s[dim] - avg;
		avgdev += d * d;
	    }
	    if (dim == 0 || avgdev > bavgdev) {
		bavgdev = avgdev;
		dsplit = dim;
		bavg = avg;
	    }
	}

	// Sort samples between left and right leaves
	for (i = 0; i < nbs; i++) {
	    if (samples[i]->s[dsplit] <= bavg) nbs_left++;
	    else break;
	}
	for (; i < nbs; i++) {
	    if (samples[i]->s[dsplit] <= bavg) {
		p = samples[nbs_left];
		samples[nbs_left++] = samples[i];
		samples[i] = p;
	    }
	}
	if (nbs_left == 0 || nbs_left == nbs) leaf = 1;
    } else leaf = 1;
    if (leaf) {
	// This is a leaf
	parent->count = nbs;
	parent->wgtCent = **samples;
	if (nbs != 1) for (i = 0; i < nb_dim; i++) parent->wgtCent.s[i] *= nbs;
	parent->point = *samples;
	parent->right = NULL;
	return parent + 1;	
    } else {
	parent->value = bavg;
	parent->dim = dsplit;
	// Call the recursion for the left part of the kdtree
	parent->right = camColorKdTreeRecurse(parent + 1, samples, nbs_left, nb_dim);
	parent->wgtCent = (parent + 1)->wgtCent;
	parent->count = (parent + 1)->count;
	// Call the recursion for the right part of the kdtree
	retval = camColorKdTreeRecurse(parent->right, samples + nbs_left, nbs - nbs_left, nb_dim);
	for (i = 0; i < nb_dim; i++) parent->wgtCent.s[i] += parent->right->wgtCent.s[i];
	parent->count += parent->right->count;
	return retval;
    }
}

CamColorKdTreeNode *camColorKdTreeBuild(CamColorSample **samples, int nb_samples, int nb_dim)
{
    CamColorKdTreeNode *kdtree = (CamColorKdTreeNode*)malloc(sizeof(CamColorKdTreeNode) * nb_samples * 2);
    camColorKdTreeRecurse(kdtree, samples, nb_samples, nb_dim); 
    return kdtree;
}

int camColorKdTreeDotRecurse(CamColorKdTreeNode *n, FILE *handle, int id)
{
    char s[64];
    int id2;
    if (IS_LEAF(n)) {
	// This is a leaf
	sprintf(s, "%d, %d (%d)", n->point->s[0], n->point->s[1], n->count);
	fprintf(handle, "\tn%d [label=\"%s\", shape=box];\n", id, s);
	return id + 1;
    } else {
	fprintf(handle, "\tn%d [label=\"s[%d] > %d\"];\n", id, n->dim, n->value); 
	id2 = camColorKdTreeDotRecurse(n + 1, handle, id + 1);
	fprintf(handle, "\tn%d->n%d [label=\"F\"];\n", id, id + 1);
	fprintf(handle, "\tn%d->n%d [label=\"T\"];\n", id, id2);
	return camColorKdTreeDotRecurse(n->right, handle, id2);
    }
}

int camColorKdTreeDot(CamColorKdTreeNode *root, const char *filename)
{
    FILE *handle;
    handle = fopen(filename, "wt");
    if (handle == NULL) return 0;
    fprintf(handle, "digraph kdtree\n{\n\trankdir=LR;\n");
    camColorKdTreeDotRecurse(root, handle, 0);
    fprintf(handle, "}\n");
    fclose(handle);
    return 1;
}

void camColorKdTreeSVGRecurse(CamColorKdTreeNode *n, FILE *handle, int *boundaries)
{
    int color;
    int cb;
    if (IS_LEAF(n)) {
	// This is a leaf
	color = 255;
	fprintf(handle, "<rect x=\"%d\" y=\"%d\" width=\"%d\" height=\"%d\" style=\"fill:rgb(%d,%d,%d);stroke-width:1;stroke:rgb(0,0,0)\"/>\n", 
		boundaries[0], boundaries[2],
		boundaries[1] - boundaries[0],
		boundaries[3] - boundaries[2],
	        color, color, color);
    } else {
	cb = boundaries[n->dim * 2 + 1];
	boundaries[n->dim * 2 + 1] = n->value;
	// Call the recursion for the left part of the kdtree
	camColorKdTreeSVGRecurse(n + 1, handle, boundaries);
	boundaries[n->dim * 2 + 1] = cb;
	cb = boundaries[n->dim * 2];
	boundaries[n->dim * 2] = n->value;
	// Call the recursion for the right part of the kdtree
	camColorKdTreeSVGRecurse(n->right, handle, boundaries);
	// We reset the boundaries so that it is kept untouched by this recursive function
	boundaries[n->dim * 2] = cb;
    }
}

void camColorKdTreeSVGRecurse2(CamColorKdTreeNode *n, FILE *handle)
{
    if (IS_LEAF(n)) {
	// This is a leaf
	fprintf(handle, "<circle cx=\"%d\" cy=\"%d\" r=\"1\" stroke=\"red\" fill=\"red\" fill-opacity=\"0.1\"/>\n",
		n->point->s[0], n->point->s[1]);
    } else {
	camColorKdTreeSVGRecurse2(n + 1, handle);
	camColorKdTreeSVGRecurse2(n->right, handle);
    }
}

int camColorKdTreeSVG(CamColorKdTreeNode *root, const char *filename)
{
    FILE *handle;
    int i, boundaries[NB_DIMENSIONS * 2];
    for (i = 0; i < NB_DIMENSIONS; i++) {
	boundaries[i * 2] = 0;
	boundaries[i * 2 + 1] = 256;
    }
    handle = fopen(filename, "wt");
    if (handle == NULL) return 0;
    fprintf(handle, "<?l version=\"1.0\" standalone=\"no\"?>\n<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
    fprintf(handle, "<svg width=\"256\" height=\"256\" version=\"1.1\" lns=\"http://www.w3.org/2000/svg\">\n");
    camColorKdTreeSVGRecurse(root, handle, boundaries);
    camColorKdTreeSVGRecurse2(root, handle);
    fprintf(handle, "</svg>\n");
    fclose(handle);
    return 1;
}

int camColorKdTreeDist(CamColorSample *s1, CamColorSample *s2, int nb_dim)
{
    int i, d, e;
    e = s1->s[0] - s2->s[0];
    d = e * e;
    for (i = 1; i < nb_dim; i++) {
	e = s1->s[i] - s2->s[i];
	d += e * e;
    }
    return d;
}
	    
int camColorKdTreeFilterIsFarther(CamColorCenter *z, CamColorCenter *zstar, int C[], int nb_dim)
{
    int i, distz, distzstar, d;
    CamColorSample vH;
    for (i = 0; i < nb_dim; i++) {
	if (z->point.s[i] > zstar->point.s[i]) {
	    vH.s[i] = C[(i << 1) + 1];
	} else {
	    vH.s[i] = C[i << 1];
	}
    }
    d = vH.s[0] - z->point.s[0];
    distz = d * d;
    d = vH.s[0] - zstar->point.s[0];
    distzstar = d * d;
    for (i = 1; i < nb_dim; i++) {
	d = vH.s[i] - z->point.s[i];
	distz += d * d;
	d = vH.s[i] - zstar->point.s[i];
	distzstar += d * d;
    }
    return (distz > distzstar);
}

void camColorKdTreeFilter(CamColorKdTreeNode *u, int C[], CamColorCenter *Z[], int sz, int nb_dim)
{
    int i, j, mind, szc, zstar;
    CamColorCenter *z, *tmp;
    CamColorSample midpoint;

    if (IS_LEAF(u)) { // if u is a leaf
	// z* <- the closest point in Z to u.point
	mind = camColorKdTreeDist(&Z[0]->point, u->point, nb_dim);
	zstar = 0;
	for (i = 1; i < sz; i++) {
	    j = camColorKdTreeDist(&Z[i]->point, u->point, nb_dim);
	    if (j < mind) {
		mind = j;
		zstar = i;
	    }
	}
	// z*.wgtCent <- zstart.wgtCent + u.point
	for (i = 0; i < nb_dim; i++)
	    Z[zstar]->wgtCent.s[i] += u->wgtCent.s[i];
	// z*.count <- z*.count + 1
	Z[zstar]->count += u->count;
    } else {
	// z* <- the closest point in Z to C's midpoint
	for (i = 0; i < nb_dim; i++)
	    midpoint.s[i] = (C[i << 1] + C[(i << 1) + 1]) >> 1;
	mind = camColorKdTreeDist(&Z[0]->point, &midpoint, nb_dim);
	zstar = 0;
	for (i = 1; i < sz; i++) {
	    j = camColorKdTreeDist(&Z[i]->point, &midpoint, nb_dim);
	    if (j < mind) {
		mind = j;
		zstar = i;
	    }
	}
	tmp = Z[0]; Z[0] = Z[zstar]; Z[zstar] = tmp; // zstar is put in front
	// for each (z belonging to Z \ {z*})
	for (szc = sz, i = 1; i < szc; i++) {
	    z = Z[i];
	    // if (z.isFarther(z*,C)) Z <- Z \ {z}
	    if (camColorKdTreeFilterIsFarther(z, Z[zstar], C, nb_dim)) {
		tmp = Z[szc - 1];
		Z[szc - 1] = z;
		Z[i] = tmp;
		szc--; i--; 
	    }
	}
	// if (|Z| = 1) {
	if (szc == 1) {
	    // z*.wgtCent <- z*.wgtCent + u.wgtCent
	    for (i = 0; i < nb_dim; i++)
		Z[zstar]->wgtCent.s[i] += u->wgtCent.s[i];
	    // z*.count <- z*.count + u.count
	    Z[zstar]->count += u->count;
	} else {
	    // Filter(u.left, Z);
	    i = C[(u->dim << 1) + 1]; 
	    C[(u->dim << 1) + 1] = u->value;
	    camColorKdTreeFilter(u + 1, C, Z, szc, nb_dim);
	    C[(u->dim << 1) + 1] = i;
	    // Filter(u.right, Z);
	    i = C[u->dim << 1]; 
	    C[u->dim << 1] = u->value + 1;
	    camColorKdTreeFilter(u->right, C, Z, szc, nb_dim);
	    C[u->dim << 1] = i;
	}
    }
}

void test1()
{
    int i;
#define TEST1_NB_SAMPLES 10 
    CamColorSample samples[TEST1_NB_SAMPLES], *pt[TEST1_NB_SAMPLES];
    CamColorKdTreeNode *kdtree;
    
    for (i = 0; i < TEST1_NB_SAMPLES; i++) {
	samples[i].s[0] = rand() % 256;
	samples[i].s[1] = rand() % 256;
	pt[i] = &samples[i];
    }
    kdtree = camColorKdTreeBuild(pt, TEST1_NB_SAMPLES, 2);
    camColorKdTreeDot(kdtree, "test1.dot");
    camColorKdTreeSVG(kdtree, "test1.svg");
    free(kdtree);
}

void test2()
{
    int i, j, k, l;
#define TEST2_NB_CENTERS 5 
#define TEST2_NB_SAMPLES 100 
#define TEST2_EXTENT 80
    CamColorSample samples[TEST2_NB_SAMPLES * TEST2_NB_CENTERS], *pt[TEST2_NB_SAMPLES * TEST2_NB_SAMPLES];
    CamColorKdTreeNode *kdtree;
    CamColorSample centers[TEST2_NB_CENTERS] = {
	{100, 50}, {50, 200}, {110, 220}, {150, 180}, {220, 70}
    };
    
    for (i = 0, k = 0; i < TEST2_NB_CENTERS; i++) {
	for (j = 0; j < TEST2_NB_SAMPLES; j++, k++) {
	    for (l = 0; l != 2; l++) {
		do samples[k].s[l] = centers[i].s[l] + (rand() % TEST2_EXTENT) - TEST2_EXTENT / 2;
		while (samples[k].s[l] < 0 || samples[k].s[l] > 255);
	    }
	    pt[k] = &samples[k];
	}
    }
    kdtree = camColorKdTreeBuild(pt, TEST2_NB_SAMPLES * TEST2_NB_CENTERS, 2);
    camColorKdTreeDot(kdtree, "test2.dot");
    camColorKdTreeSVG(kdtree, "test2.svg");
    free(kdtree);
}

int main()
{
    test1();
    test2();
    return 0;
}

    
