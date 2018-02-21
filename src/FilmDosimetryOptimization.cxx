#include "FilmDosimetryOptimization.h"

/*!
    \class FilmDosimetryOptimization
    \brief The FilmDosimetryOptimization class encapsulates
    robust multi-channel optimization in EBT2/3 dosimetry.
    \since 0.


    The aim of this class is to implement the multichannel method.

    The metodology is described at:
    \sa
   https://www.sciencedirect.com/science/article/pii/S0010465513000805?via%3Dihub
*/

/*!
    \fn FilmDosimetryOptimization::FilmDosimetryOptimization(cv::String
   &fileName, std::vector<int> eqt, std::vector<double> pol)

    Constructs a optimization object from \a tiff file path, \a array of
   qualibration equation index and \a polynomial coefficients for RGB channels
*/
FilmDosimetryOptimization::FilmDosimetryOptimization(cv::String fileName,
                                                     std::vector<int> eqt,
                                                     std::vector<double> pol)
    : m_fileName(fileName), m_eqt(eqt), m_pol(pol) {

  readImageTiff16bits(m_fileName);
  int nRows = m_image.size[0];
  int nCols = m_image.size[1];

  // create empty m_Emap matrix
  m_Emap = cv::Mat(nRows, nCols, CV_64FC1);

  // create empty dose matrix RGB
  m_DoseMap = cv::Mat(nRows, nCols, CV_64FC3);
}

FilmDosimetryOptimization::~FilmDosimetryOptimization() {}

/*!
    \fn FilmDosimetryOptimization::readImageTiff16bits(cv::String &fileName)

    Reads the 48bits RBG tiff file from \a file path *.tiff
    It is called by the constructor to set member variables

*/
void FilmDosimetryOptimization::readImageTiff16bits(cv::String &fileName) {

  auto image = cv::imread(fileName, cv::IMREAD_ANYDEPTH | cv::IMREAD_ANYCOLOR);

  // OpenCV reads using BGR channel order
  auto fm = cv::Mat(image.size(), CV_64FC3);
  image.convertTo(fm, CV_64FC3);
  // converting 16bit opencv Matrix to OD
  double normFactor = 65535.0;
  int nRows = image.size[0];
  int nCols = image.size[1];
  for (int y = 0; y < nRows; y++) {
    for (int x = 0; x < nCols; x++) {
      // OpenCV reads using BGR channel order
      auto intensity = fm.at<cv::Vec3d>(y, x);
      auto blue = intensity.val[0] / normFactor;
      auto green = intensity.val[1] / normFactor;
      auto red = intensity.val[2] / normFactor;

      // swap channels to RGB order
      intensity.val[0] = -log10(red);
      intensity.val[1] = -log10(green);
      intensity.val[2] = -log10(blue);

      fm.at<cv::Vec3d>(y, x) = intensity;
    }
  }
  // copy to member variables
  image.copyTo(m_image);
  fm.copyTo(m_image_od);
}

/*!
    \fn double FilmDosimetryOptimization::inverseEquation(double x, int eqt,
                                                  std::vector<double> p)
    Returns \c Dose (cGy) for a given \a optical density, \a equation index and
   \a polynomial coefficientes

   It is inverse equation because the polynomial coefficients were generated
   using:
    x-axis : Dose in cGy
    y-axis : OD

   eq 1 : OD = a3 * log(D/A + 1)^2 + a2 * log(D/A + 1)+ a1
   eq 2 : OD = a3 * D^2 + a2 * D + a1
   eq 3 : OD = a3 * atan(D/B)^2 + a2 * atan(D/B) + a1

   Where:
     A = 100.0
     B = 500.0

   D = Dose in cGy
   OD = Optical density

*/
double
FilmDosimetryOptimization::inverseEquation(const double &od, const int &eqt,
                                           const std::vector<double> &p) {

  double A = 100.0;
  double B = 500.0;
  double Delta = 4 * p[2] * od - 4 * p[0] * p[2] + pow(p[1], 2);
  // check real root
  if (od > Delta && Delta > 0.) {
    if (eqt == 1) {
      auto param = sqrt(Delta) / (2 * p[2]) - p[1] / (2 * p[2]);
      return (exp(param) - 1) * A;
    }
    if (eqt == 2) {
      return (sqrt(Delta) - p[1]) / (2 * p[2]);
    }
    if (eqt == 3) {
      auto frac = (sqrt(Delta) - p[1]) / (2 * p[2]);
      return tan(frac) * B;
    }
  } else {
    return 0.;
  }
  return 0.;
}

/*!
    \fn double FilmDosimetryOptimization::diffCalibrationEquation(double x, int
   eqt, std::vector<double> p)


   Returns \c derivative of calibration equations for a given \a optical
   density, \a equation index and \a polynomial coefficientes

   Calibration equations

   eq 1 : OD = a3 * log(D/A + 1)^2 + a2 * log(D/A + 1)+ a1
   eq 2 : OD = a3 * D^2 + a2 * D + a1
   eq 3 : OD = a3 * atan(D/B)^2 + a2 * atan(D/B) + a1

   Where:
     A = 100.0
     B = 500.0

   D = Dose in cGy
   OD = Optical density

*/
double FilmDosimetryOptimization::diffCalibrationEquation(
    const double &dose, const int &eqt, const std::vector<double> &p) {
  double A = 100.0;
  double B = 500.0;

  if (eqt == 1) {
    return (2 * p[2] * log((A + dose) / A) + p[1]) / (A + dose);
  }
  if (eqt == 2) {
    return 2 * p[2] * dose + p[1];
  }
  if (eqt == 3) {
    return ((2 * p[2] * atan(dose / B) + p[1]) * B) /
           (pow(B, 2) + pow(dose, 2));
  }
  return 0.;
}

/*!
    \fn double FilmDosimetryOptimization::objectiveFunction(double dz,
   cv::Vec3d od, std::vector<int> eqt, std::vector<double> pol


   Helper Derivative of Objective function to find equation roots
   Function minimum is where derivative equals 0.;

   Returns Diff of objective function.

   Receives \a dimensionless factor dz, \a od vector containing pixels values
   RGB \a equation vector RBG, \a contiguous array containing all calibration
   coefficients RGB (size 9)

*/
double FilmDosimetryOptimization::objectiveFunction(
    const double &dz, const cv::Vec3d &od, const std::vector<int> &eqt,
    const std::vector<double> &pr, const std::vector<double> &pg,
    const std::vector<double> &pb) {
  /*
   * Helper Derivative of Objective function to find equation roots
   * Function minimum is where derivative equals 0.;
   */

  // calculating "disturbed" od
  double odr = od[0] * dz;
  double odg = od[1] * dz;
  double odb = od[2] * dz;

  // getting polinomial coefficients
  //  std::vector<double> pr = {pol[0], pol[1], pol[2]};
  //  std::vector<double> pg = {pol[3], pol[4], pol[5]};
  //  std::vector<double> pb = {pol[6], pol[7], pol[8]};

  // getting doses
  auto Dr = inverseEquation(odr, eqt[0], pr);
  auto Dg = inverseEquation(odg, eqt[1], pg);
  auto Db = inverseEquation(odb, eqt[2], pb);

  auto dODr = diffCalibrationEquation(Dr, eqt[0], pr);
  auto dODg = diffCalibrationEquation(Dg, eqt[1], pg);
  auto dODb = diffCalibrationEquation(Db, eqt[2], pb);

  // inverse function derivative
  dODr = 1 / dODr;
  dODg = 1 / dODg;
  dODb = 1 / dODb;

  //# Chain rule derivative
  auto dDr_da = dODr * od[0];
  auto dDg_da = dODg * od[1];
  auto dDb_da = dODb * od[2];

  // differential of Robust objective function

  auto sc = 1 / 9.0;

  auto Drg = Dr - Dg;
  auto Dgb = Dg - Db;
  auto Dbr = Db - Dr;

  auto omega = (Drg / sqrt(pow(Drg, 2) + sc)) * (dDr_da - dDg_da) / 3.0 +
               (Dgb / sqrt(pow(Dgb, 2) + sc)) * (dDg_da - dDb_da) / 3.0 +
               (Dbr / sqrt(pow(Dbr, 2) + sc)) * (dDb_da - dDr_da) / 3.0;

  return omega;
}

/*!
    \fn void FilmDosimetryOptimization::optimizeDosimetry()

   Zeroin equation solver - finds zeros of nonlinear equations with no
   derivatives requirement of Jacobians like Newton's solver.

    Adapted from:

        ref. htps://en.wikipedia.org/wiki/Brent%27s_method

*/
void FilmDosimetryOptimization::optimizeDosimetry(
    const std::vector<double> pr, const std::vector<double> pg,
    const std::vector<double> pb) {
  /*
      Zeroin equation solver - finds zeros of nonlinear equations with no
     derivatives

      Adapted from:

          ref. htps://en.wikipedia.org/wiki/Brent%27s_method

  */
  double tol = 2.2204460492503131e-16;
  // getting polinomial coefficients
  int nRows = m_image.size[0];
  int nCols = m_image.size[1];
  for (int i = 0; i < nRows; i++) {   // y
    for (int j = 0; j < nCols; j++) { // x

      double x = 1.0;
      auto intensity = m_image_od.at<cv::Vec3d>(i, j);
      auto fx = objectiveFunction(x, intensity, m_eqt, pr, pg, pb);
      if (cmpf(fx, 0.)) {
        m_Emap.at<double>(i, j) = x;
      } else {
        double dx = 1 / 50.0;
        //# Find change of  sign.
        double twosqrt = sqrt(2);
        double a = x;
        double fa = fx;
        double b = x;
        double fb = fx;

        while (fa > 0. == fb > 0.) {
          /* # for isig in xrange(loop):
          # if not ((fa > 0) == (fb > 0)):
          # break*/
          dx *= twosqrt;
          a = x - dx;
          fa = objectiveFunction(a, intensity, m_eqt, pr, pg, pb);

          if (fa > 0. != fb > 0.) { // check for different sign
            break;
          }
          b = x + dx;
          fb = objectiveFunction(b, intensity, m_eqt, pr, pg, pb);
        }

        auto fc = fb;
        auto c = b;
        auto d = b - a;
        auto e = d;

        while (!cmpf(fb, 0.) && !cmpf(a, b)) {

          if (fb > 0. == fc > 0) {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
          }

          auto tfc = std::fabs(fc);
          auto tfb = std::fabs(fb);

          if (tfc < tfb) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
          }

          // Convergence test and possible exit
          auto m = 0.5 * (c - b);
          auto tbb = std::fabs(b);
          auto toler = 2.0 * tol * std::max(tbb, 1.0);

          auto tmm = std::fabs(m);
          if (tmm <= toler || cmpf(fb, 0.)) {
            break;
          }
          // Choose bisection or interpolation
          auto tt1 = std::fabs(e) < toler;
          auto tt2 = std::fabs(fa) <= std::fabs(fb);
          if (tt1 || tt2) {
            //# Bisection
            d = m;
            e = m;
          } else {
            //#  Interpolation
            //                        double p = 0.0, q = 0.0;
            double p, q;
            auto s = fb / fa;
            if (cmpf(a, c)) {
              // Linear interpolation
              p = 2.0 * m * s;
              q = 1.0 - s;

            } else {
              // Inverse quadratic interpolation
              q = fa / fc;
              auto r = fb / fc;
              p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0));
              q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (p > 0.) {
              q = -q;
            } else {
              p = -p;
            }

            // Is interpolated point acceptable ?
            auto t1 = 2.0 * p < (3.0 * m * q - std::fabs(toler * q));
            auto t2 = p < std::fabs(0.5 * e * q);

            if (t1 && t2) {
              e = d;
              d = p / q;
            } else {
              d = m;
              e = m;
            }
          }
          //# % Next point
          a = b;
          fa = fb;
          auto td = std::fabs(d);

          if (td > toler) {
            b = b + d;
          } else if (b > c) {
            b = b - toler;
          } else {
            b = b + toler;
          }
          fb = objectiveFunction(b, intensity, m_eqt, pr, pg, pb);
        }
        m_Emap.at<double>(i, j) = b;
      }
    }
  }
}

/*!
    \fn FilmDosimetryOptimization::calcOptimizedDoseRGB()

    Calculates the optimized film dose for all color channels using robust
    multi-channel optimization.

*/
void FilmDosimetryOptimization::calcOptimizedDoseRGB() {

  // getting polinomial coefficients
  std::vector<double> pr = {m_pol[0], m_pol[1], m_pol[2]};
  std::vector<double> pg = {m_pol[3], m_pol[4], m_pol[5]};
  std::vector<double> pb = {m_pol[6], m_pol[7], m_pol[8]};

  // calculate disturbance map
  optimizeDosimetry(pr, pg, pb);

  // calculate optimized dose pixel by pixel
  int nRows = m_image_od.size[0];
  int nCols = m_image_od.size[1];

  for (int y = 0; y < nRows; y++) {
    for (int x = 0; x < nCols; x++) {

      auto intensity = m_image_od.at<cv::Vec3d>(y, x);
      auto disturbance = m_Emap.at<double>(y, x);

      // corrected od per channel
      auto red = intensity.val[0] * disturbance;
      auto green = intensity.val[1] * disturbance;
      auto blue = intensity.val[2] * disturbance;

      // calculating dose per channel (using inverse calibration
      // equation equations)
      auto Dr = inverseEquation(red, m_eqt[0], pr);
      auto Dg = inverseEquation(green, m_eqt[1], pg);
      auto Db = inverseEquation(blue, m_eqt[2], pb);

      cv::Vec3d doseRGB;
      doseRGB.val[0] = Dr;
      doseRGB.val[1] = Dg;
      doseRGB.val[2] = Db;

      m_DoseMap.at<cv::Vec3d>(y, x) = doseRGB;
    }
  }
}
/*!
    \fn  FilmDosimetryOptimization::printResultMatrix()

    Print disturbance map and Dose for all channels to terminal.
*/
void FilmDosimetryOptimization::printResultMatrix() {
  std::cout << "Disturbance map: " << std::endl;
  std::cout << cv::format(m_Emap, cv::Formatter::FMT_NUMPY) << std::endl;
  std::cout << "Dose map cGy: " << std::endl;
  std::cout << cv::format(m_DoseMap, cv::Formatter::FMT_NUMPY) << std::endl;
}

cv::Mat FilmDosimetryOptimization::getImage() const { return m_image; }

cv::Mat FilmDosimetryOptimization::getDoseMap() const { return m_DoseMap; }

/*!
    \property FilmDosimetryOptimization::Emap()
    \brief Disturbance map after optimization
 */
cv::Mat FilmDosimetryOptimization::Emap() const { return m_Emap; }

/*!
    \property FilmDosimetryOptimization::Emap()
    \brief Disturbance map after optimization
 */
cv::Mat FilmDosimetryOptimization::image_od() const { return m_image_od; }
