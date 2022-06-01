#ifndef __MAX_SEARCH_H__
#define __MAX_SEARCH_H__

template<int N>
struct max_s {
  template<typename T>
  static T max(T *a) {
    T m0 = max_s<N/2>::max(a);
    T m1 = max_s<N-N/2>::max(a+N/2);

    return m0 > m1 ? m0 : m1;
  }
};

template<> 
struct max_s<1> {
  template<typename T>
  static T max(T *a) {
    return a[0];
  }
};

template<int N, typename T>
T max(T *a) {
  return max_s<N>::max(a);
}

#endif

