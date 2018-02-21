#ifndef FILMDOSIMETRYOPTIMIZATION_H
#define FILMDOSIMETRYOPTIMIZATION_H

#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <opencv2/opencv.hpp>
#include <vector>

class FilmDosimetryOptimization {
public:
  FilmDosimetryOptimization(cv::String fileName, std::vector<int> eqt,
                            std::vector<double> pol);
  virtual ~FilmDosimetryOptimization();

private:
  // filename
  cv::String m_fileName;

  // original image 16 bits BGR
  cv::Mat m_image;
  // OD image RGB
  cv::Mat m_image_od;
  // Disturbance map
  cv::Mat m_Emap;
  // dose map (RGB)
  cv::Mat m_DoseMap;

  // calibration equation vector RGB
  std::vector<int> m_eqt;
  // calibration equation vector RGB
  std::vector<double> m_pol;

public:
  /*
    Calculate dose using robust MC - optimization
  */
  void calcOptimizedDoseRGB();
  void printResultMatrix();
  // Getters
  cv::Mat image_od() const;
  cv::Mat Emap() const;
  cv::Mat getImage() const;
  cv::Mat getDoseMap() const;

private:
  // reads a tiff 48 bits file to optical density matrix
  void readImageTiff16bits(cv::String &fileName);
  /*!
      Inverse funtion, return Dose (cGy) per Optical Density.
  :param x: Optical Density OD
  :param eqt:  inverse equation type:
          1 : 'y = Sigma a_nln(x/A + 1)^(^n^-^1^)'
          2 : 'y = Sigma a_nx^(^n^-^1^)'
          3 : 'y = Sigma a_natan(x/B)^(^n^-^1^)'

  :param p: polinomial coeficients ( curve fit ) p=flipud(p);
  :return: Dose ( cGy )
  */
  double inverseEquation(const double &od, const int &eqt,
                         const std::vector<double> &p);

  double diffCalibrationEquation(const double &dose, const int &eqt,
                                 const std::vector<double> &p);
  /*
      Calculate the first derivative of calibration curve equations.

      x-axis : Dose in cGy
      y-axis : Optical Density

  :param x: Calibration doses (cGy) ( Array or scalar)
  :param eqt:  equation type:
          1 : 'y = \Sigma a_nln(x/A + 1)^(^n^-^1^)'
          2 : 'y = \Sigma a_nx^(^n^-^1^)'
          3 : 'y = \Sigma a_natan(x/B)^(^n^-^1^)'
  :param p: polinomial coeficients ( curve fit ) p=flipud(p);
  :return: First derivative
  */
  double objectiveFunction(const double &dz, const cv::Vec3d &od,
                           const std::vector<int> &eqt,
                           const std::vector<double> &pr,
                           const std::vector<double> &pg,
                           const std::vector<double> &pb);
  /*
   Derivative from Robust objective function on robust multi-channel
  optimization
      ref: http://adsabs.harvard.edu/abs/2013CoPhC.184.1708A

  :param dz: optimization parameter
  :param od: optical density
  :param eqt: calibration equation index RGB (0, 1, 2)
  :param pol: polynomial fit coefficients
  vector indices:
          (0 1 2) Red
          (3 4 5) Green
          (6 7 8) Green
  :return: function eval
  */

  void optimizeDosimetry();
  /*
      Zeroin equation solver - finds zeros of nonlinear equations with no

      derivatives ref. htps://en.wikipedia.org/wiki/Brent%27s_method

  :param od: Image in OD (channel: RGB)
  :param eqt: calibration equation index
  :param pol: polynomial fit coefficients
  :return: Root of equation

  */
  /*
      Calculate dose using robust MC - optimization
  */
  bool cmpf(const double &A, const double &B,
            const double &epsilon = 2.2204460492503131e-16) {
    /* Helper method to compare double point precision numbers */
    return (std::fabs(A - B) < epsilon);
  }
};

#endif // FILMDOSIMETRYOPTIMIZATION_H
