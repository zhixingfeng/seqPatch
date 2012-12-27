#ifndef MATHOPT_H
#define MATHOPT_H

double max_elt(double *data, int len)
{
	double max_value = data[0];
	for (int i=1;i<len;i++){
		if (data[i]>max_value)
			max_value = data[i];
	}
	return max_value;
}
int max_elt(int *data, int len)
{
        int max_value = data[0];
        for (int i=1;i<len;i++){
                if (data[i]>max_value)
                        max_value = data[i];
        }
        return max_value;
}

double min_elt(double *data, double len)
{
        double min_value = data[0];
        for (int i=1;i<len;i++){
                if (data[i]<min_value)
                        min_value = data[i];
        }
        return min_value;
}

int min_elt(int *data, int len)
{
        int min_value = data[0];
        for (int i=1;i<len;i++){
                if (data[i]<min_value)
                        min_value = data[i];
        }
        return min_value;
}



#endif




