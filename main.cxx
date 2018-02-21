#include "FilmDosimetryOptimization.h"
#include <math.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>

void tic(double &t) { t = (double)cv::getTickCount(); };

double toc(double &t) {
  return ((double)cv::getTickCount() - t) / cv::getTickFrequency();
};

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("usage: Film2DoseOptimizaion.exe <Image_Path>\n");
    return -1;
  }

  // sample calibration equations
  // calibration equations RBG
  std::vector<int> eqt{1, 1, 1};

  // as contiguous array
  std::vector<double> poly{
      0.2594476,  0.15756712, 0.02068541,  // Red calibration coeffs
      0.28112982, 0.07294514, 0.04343196,  // Green calibration coeffs
      0.53939429, 0.0152481,  0.03212944}; // Blue calibration coeffs

  auto fileName = static_cast<cv::String>(argv[1]);
  auto *optimizer = new FilmDosimetryOptimization(fileName, eqt, poly);

  // run optimization
  double t;
  tic(t);
  optimizer->calcOptimizedDoseRGB();
  std::cout << "Elapsed (s): " << toc(t) << std::endl;

  //  optimizer->printResultMatrix();

  auto doseMap = optimizer->getDoseMap();
  cv::Mat rgb[3];          // destination array
  cv::split(doseMap, rgb); // split source
  cv::imshow("Dose Map Image", rgb[2]);
  cv::imshow("Disturance map Image", optimizer->Emap());
  cv::waitKey(0);

  delete optimizer;
  return 0;
}
