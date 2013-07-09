#ifndef MATH_UTILS_H
///Define this macro to prevent from including this header file more than once.
#define MATH_UTILS_H

#include "stl.h"
#include "type.h"

const double epsilon = 1e-10;

inline int round_double(double x){
	if (x > 0) return int(x + 0.5);
	else return int(x-0.5);
}

inline bool is_int(double d) {
	return (fabs(d - round_double(d)) < epsilon);
}

inline double rand_double(){
	return ((double)(rand()))/RAND_MAX;
}

class interval{
public:
	double start, end;
	interval(const interval &temp) {
		start = temp.start;
		end = temp.end;
	}
	interval(const double start, const double end){
		this->start = start;
		this->end = end;
	}
	inline double length() const {
		return end - start;
	}
	inline bool is_empty() const{
		return (length() < 0);
	}
	inline void outer_union_with(const interval &temp) { //not real union, which may generate a interval_set
		if (temp.is_empty()) return;
		if (is_empty()) {
			start = temp.start;
			end = temp.end;
		}
		start = min(start, temp.start);
		end = max(end, temp.end);
	}
	inline void intersect_with(const interval &temp){
		if (is_empty()) return;
		start = max(start, temp.start);
		end = min(end, temp.end);
	}
};

class interval_set{
public:
	vector<double> end_points;
	interval_set(){
		end_points.clear();
	}

	interval_set(const interval_set &temp) {
		end_points = temp.end_points;
	}

	interval_set(const interval &temp) {
		end_points.clear();
		end_points.push_back(temp.start);
		end_points.push_back(temp.end);
	}

	interval_set(const double start, const double end) {
		end_points.clear();
		end_points.push_back(start);
		end_points.push_back(end);
	}

	interval_set(const vector<pair<int, int> > pairs) {
		for (int i = 0; i < (int)pairs.size(); i++) {
			end_points.push_back(pairs[i].first);
			end_points.push_back(pairs[i].second);
		}
	}

	inline bool convert_to_int_pairs(vector<pair<int, int> > &pairs) {
		if (!check_int()) return false;
		for (int i = 0; i < (int)end_points.size(); i += 2) {
			pair<int, int> temp_pair(round_double(end_points[i]), round_double(end_points[i+1]));
			pairs.push_back(temp_pair);
		}
		return true;
	}

	inline bool check_valid() {
		if ((int)end_points.size() %2 != 0) return false;
		for (int i = 1; i < (int)end_points.size(); i++) {
			if (end_points[i] + epsilon < end_points[i-1]) return false;
		}
		return true;
	}

	inline bool check_int() {
		if (!check_valid()) return false;
		for (int i = 0; i < (int)end_points.size(); i++) {
			if (!is_int(end_points[i])) return false;
		}
		return true;
	}

	inline double length() {
		double total = 0;
		for (int i = 0; i < (int)end_points.size(); i += 2) {
			total += end_points[i+1] - end_points[i];
		}
		return total;
	}

	inline interval bracket() const{
		if (end_points.size() > 1) return interval(end_points[0], end_points[end_points.size() - 1]);
		else return interval(0, -1);
	}

	inline bool has_empty_intervals() {
		for (int i = (int)end_points.size() - 1; i > 0; i -= 2) {
			if (end_points[i] <= end_points[i - 1] + epsilon) {
				return true;
			}
		}
		return false;
	}

	inline void remove_empty_intervals() {
		for (int i = (int)end_points.size() - 1; i > 0; i -= 2) {
			if (end_points[i] <= end_points[i - 1] + epsilon) {
				end_points.erase(end_points.begin() + i - 1, end_points.begin() + i + 1);
			}
		}
	}

	inline void remove_redundant_inner_points(){
		for (int i = (int)end_points.size() - 1; i > 0; i--) {
			if (end_points[i] <= end_points[i - 1] + epsilon) {
				end_points.erase(end_points.begin() + i - 1, end_points.begin() + i + 1);
				i--;
			}
		}
	}

	inline void complement_with(const interval &temp) {
		interval bracket1 = bracket();
		if (end_points.size() > 0 && (temp.start > bracket1.start || temp.end < bracket1.end)) return;
		end_points.insert(end_points.begin(), temp.start);
		end_points.push_back(temp.end);
		remove_redundant_inner_points();
	}

	inline void intersect_with(const interval_set &temp) {
		interval bracket1 = bracket();
		bracket1.outer_union_with(temp.bracket());
		if (bracket1.is_empty()) return;
		complement_with(bracket1);
		interval_set temp1 = temp;
		temp1.complement_with(bracket1);
		union_with(temp1);
		complement_with(bracket1);
	}

	inline void union_with(const interval_set &temp);

	inline void break_union_with(const interval &temp) {
		int i = 0;
		for (i = 0; i < (int)end_points.size(); i++) {
			if (temp.start < end_points[i] - epsilon) break;
		}
		bool overlap = (i > 0 && fabs(temp.start - end_points[i-1]) < epsilon);
		bool inside = (i % 2 == 1);
		if (!inside) {
			end_points.insert(end_points.begin() + i, temp.start);
			i += 1;
		} else if (!overlap) {
			end_points.insert(end_points.begin() + i, 2, temp.start);
			i += 2;
		}
		while (true) {
			if (i == (int)end_points.size()) break;
			if (temp.end <= end_points[i] + epsilon) break;
			double temp_pos = end_points[i];
			if (i == (int)end_points.size() - 1 || fabs(end_points[i+1] - end_points[i]) > epsilon) end_points.insert(end_points.begin() + i, temp_pos);
			i += 2;
		}
		overlap = (i < (int)end_points.size() && fabs(temp.end - end_points[i]) < epsilon);
		inside = (((int)end_points.size() - i) % 2 == 1);
		if (!inside) {
			end_points.insert(end_points.begin() + i, temp.end);
			i += 1;
		} else if (!overlap) {
			end_points.insert(end_points.begin() + i, 2, temp.end);
			i += 2;
		}
		__ASSERT(!has_empty_intervals(), "internal error: has_empty_intervals().\n");
	}

	inline void break_union_with(const interval_set &temp) {
		for (int i = 0; i < (int)temp.end_points.size(); i += 2) {
			break_union_with(interval(temp.end_points[i], temp.end_points[i+1]));
		}
	}
};

inline void interval_set::union_with(const interval_set &temp){
	vector<pair<double, bool> > points;
	int i;
	for (i = 0; i < (int)end_points.size(); i++) {
		points.push_back(pair<double, bool>(end_points[i], i%2==0));
	}
	for (i = 0; i < (int)temp.end_points.size(); i++) {
		points.push_back(pair<double, bool>(temp.end_points[i], i%2==0));
	}
	sort(points.begin(), points.end());
	end_points.clear();
	int layer = 0;
	for (i = 0; i < (int)points.size(); i++) {
		if (points[i].second) {
			layer++;
			if (layer == 1) end_points.push_back(points[i].first);
		} else {
			layer--;
			if (layer == 0) end_points.push_back(points[i].first);
		}
	}
	remove_redundant_inner_points();
}

class folding{
private:
	double equal_gap_length;
	double gap_coef, interval_coef;
public:
	vector<double> points;
	vector<double> mapped_points;
	interval_set intervals;
	bool begin_with_gap;
	double start, end;
	double mapped_start, mapped_end;
	bool equal_gap; //false: ratio_gap
	double gap_ratio;

	inline folding() {
	}
	
	inline folding(const interval_set &intervals, const double start, const double end, const double mapped_start, const double mapped_end, const bool equal_gap = true, const double gap_ratio = 0.1){
		configure(intervals, start, end, mapped_start, mapped_end, equal_gap, gap_ratio);
	}

	inline void configure(const interval_set &intervals, const double start, const double end, const double mapped_start, const double mapped_end, const bool equal_gap = true, const double gap_ratio = 0.1){
		if (start >= end - epsilon) panic("internal error: start > end.");

		this->intervals = intervals;
		this->start = start;
		this->end = end;
		this->mapped_start = mapped_start;
		this->mapped_end = mapped_end;
		this->equal_gap = equal_gap;
		this->gap_ratio = gap_ratio;
		
		this->intervals.intersect_with(interval_set(start, end));
		points = this->intervals.end_points;
		begin_with_gap = false;
		if (points.size() == 0 || start < points[0] - epsilon) {
			begin_with_gap = true;
			points.insert(points.begin(), start);
		}
		if (end > points[points.size() - 1] + epsilon) {
			points.push_back(end);
		}

		double total_interval = 0;
		double total_gap = 0;
		for (int i = 0; i < (int)points.size() - 1; i++) {
			if (begin_with_gap == (i%2==0)) total_gap += points[i+1] - points[i];
			else total_interval += points[i+1] - points[i];
		}

		double mapped_length = mapped_end - mapped_start;
		double mapped_interval_length = mapped_length * (1 - gap_ratio);
		double mapped_gap_length = mapped_length - mapped_interval_length;
		int num_gaps = (int)points.size()/2;
		if (!begin_with_gap && (points.size()%2==0)) num_gaps--;

		if (num_gaps == 0) {
			mapped_gap_length = 0;
			mapped_interval_length = mapped_length;
		}

		if (points.size() == 2 && num_gaps == 1) {
			mapped_gap_length = mapped_length;
			mapped_interval_length = 0;
		}

		if (num_gaps == 0) {
			equal_gap_length = 0;
			gap_coef = 0;
		} else {
			equal_gap_length = mapped_gap_length / num_gaps;
			gap_coef = mapped_gap_length / total_gap;
		}

		if (points.size() == 2 && num_gaps == 1) {
			interval_coef = 0;
		} else {
			interval_coef = mapped_interval_length / total_interval;
		}

		mapped_points.resize(points.size());
		for (int i = 0; i < (int)points.size(); i++) {
			if (i == 0) { 
				mapped_points[i] = mapped_start;
			} else if (begin_with_gap == (i%2==1)) {
				if (equal_gap) {
					mapped_points[i] = mapped_points[i-1] + equal_gap_length;
				} else {
					mapped_points[i] = mapped_points[i-1] + (points[i] - points[i-1]) * gap_coef;
				}
			} else {
				mapped_points[i] = mapped_points[i-1] + (points[i] - points[i-1]) * interval_coef;
			}
		}
		if (fabs(mapped_points[mapped_points.size() - 1] - mapped_end) > epsilon) panic("internal error, mapping length inconsistent.");
	}

	double map(double point) {
       int low = 0, high = (int)points.size();
       while (low < high) {
           int mid = (low + high)/2;
           if (points[mid] < point)
               low = mid + 1;
           else
               high = mid; 
       }
	   if (low == (int)points.size()) {
		   if (point < end) panic("internal error, bad mapping.");
		   return mapped_end + (point - end) * interval_coef;
	   } else if (low == 0) {
		   if (point > start) panic("internal error, bad mapping.");
		   return mapped_start - (start - point) * interval_coef;
	   } else {
		   return ((point - points[low - 1]) * mapped_points[low] + (points[low] - point) * mapped_points[low - 1]) / (points[low] - points[low - 1]);
	   }
	}

	double inverse_map(double mapped_point) {
       int low = 0, high = (int)mapped_points.size();
       while (low < high) {
           int mid = (low + high)/2;
           if (mapped_points[mid] < mapped_point)
               low = mid + 1;
           else
               high = mid; 
       }
	   if (low == (int)mapped_points.size()) {
		   if (mapped_point < mapped_end) panic("internal error, bad mapping.");
		   return end;
	   } else if (low == 0) {
		   if (mapped_point > mapped_start) panic("internal error, bad mapping.");
		   return start;
	   } else {
		   return ((mapped_point - mapped_points[low - 1]) * points[low] + (mapped_points[low] - mapped_point) * points[low - 1]) / (mapped_points[low] - mapped_points[low - 1]);
	   }
	}
};

inline double sum(const vector<double> &data) {
	double result = 0;
	for (int i = 0; i < (int)data.size(); i++) {
		result += data[i];
	}
	return result;
}

inline double mean(const vector<double> &data) {
	return sum(data) / (int)data.size();
}

inline double var(const vector<double> &data) {
	double sum1 = 0, sum2 = 0;
	for (int i = 0; i < (int)data.size(); i++) {
		sum1 += data[i];
		sum2 += data[i] * data[i];
	}
	return (sum2 - sum1 * sum1 / (int)data.size()) / ((int)data.size() - 1);
}

inline double mean_c(double * data, int n)
{
	double rl = 0;
	for (int i=0;i<n;i++){
		rl += data[i];
	}
	return rl/n;
}

inline double var_c(double *data, int n)
{
	double sum1 = 0, sum2 = 0;
        for (int i = 0; i < n; i++) {
                sum1 += data[i];
                sum2 += data[i] * data[i];
        }

        return (sum2 - sum1 * sum1 / n ) / (n - 1);

}

inline double stddev(const vector<double> &data) {
	return sqrt(var(data));
}


//  The following code is written by Zhixing Feng
inline vector<double> filter_outlier(vector<double> &data, double k = 2.5, string tail="right")
{
	sort(data.begin(), data.end());
	int bd = int( data.size() * 0.25);
	double low = data[bd]; double high = data[data.size() -1 -bd];
	double len = high - low;
	vector<double> data_ft;
	if (tail=="right")
		for (int i=0;i<(int)data.size();i++) if (data[i]<=high+k*len) data_ft.push_back(data[i]);
	if (tail=="left")
                for (int i=0;i<(int)data.size();i++) if (data[i]>=low-k*len) data_ft.push_back(data[i]);
	if (tail=="both")
		for (int i=0;i<(int)data.size();i++) if (data[i]>=low-k*len && data[i]<=high+k*len) data_ft.push_back(data[i]);
	//Rprintf("%lf\n%lf\n",low,high);
	return data_ft;
}


#endif // MATH_UTILS_H
