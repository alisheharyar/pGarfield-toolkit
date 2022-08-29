#ifndef PARABOL_H
#define PARABOL_H

/*
Copyright (c) 2001 I. B. Smirnov

Permission to use, copy, modify, distribute and sell this file
and its documentation for any purpose is hereby granted without fee,
provided that the above copyright notice, this permission notice,
and notices about any modifications of the original text
appear in all copies and in supporting documentation.
It is provided "as is" without express or implied warranty.
*/

namespace Heed {

/// Solution of a quadratic equation.

class Parabola {
 public:
  double a() const { return da; }
  double b() const { return db; }
  double c() const { return dc; }
  void put_a(const double fa) {
    da = fa;
    s_det = 0;
    s_dxzero = 0;
  }
  void put_b(const double fb) {
    db = fb;
    s_det = 0;
    s_dxzero = 0;
  }
  void put_c(const double fc) {
    dc = fc;
    s_det = 0;
    s_dxzero = 0;
  }

  /// Default constructor.
  Parabola() = default;
  /// Constructor from coefficients.
  Parabola(double fa, double fb, double fc)
      : da(fa), db(fb), dc(fc) {}
  /// Constructor from three points.
  Parabola(double x[3], double y[3]);
  /// Constructor from three points. 
  /// At the third one, the derivative of the function is supplied instead of
  /// the function. 
  Parabola(double x[3], double y[3], int);
  /// Constructor from three points.
  Parabola(double x1, double x2, double x3, double y1, double y2, double y3);

  /// Copy constructor.
  Parabola(const Parabola& f);
  /// Copy assignment operator.
  Parabola& operator=(const Parabola& p) = default;

  /// Evaluate the function.
  double eval(const double x) const { return da * x * x + db * x + dc; }

  // Returns number of solutions. First is the least.
  int find_zero(double xzero[2]) const;
  double find_maxmin();

  double determinant() const {
    const Parabola& t = (*this);
    if (s_det == 0) {
      t.s_det = 1;
      t.det = db * db - 4 * da * dc;
    }
    return det;
  }

 private:
  double da = 0., db = 0., dc = 0.;
  mutable int s_det = 0;
  mutable double det = 0.;
  mutable int s_dxzero = 0;
  mutable int qdxzero = 0;
  mutable double dxzero[2];
};

std::ostream& operator<<(std::ostream& file, const Parabola& f);
}

#endif
