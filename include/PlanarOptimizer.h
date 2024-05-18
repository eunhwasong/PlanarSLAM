#ifndef PLANAR_SLAM_OPTIMIZER_H
#define PLANAR_SLAM_OPTIMIZER_H
#pragma once

#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <atomic>
#include <LoopClosingTypes.h>
#include "g2o/core/block_solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/types/sba/types_six_dof_expmap.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/types/sim3/types_seven_dof_expmap.h"



namespace EdgeSLAM {
	class MapPoint;
	class Frame;
	class KeyFrame;
	class Map;
	class SLAM;

	class PlanarOptimizer {
	public:
		void static BundleAdjustment(SLAM* SLAM, const std::vector<KeyFrame*> &vpKF, const std::vector<MapPoint*> &vpMP,
			int nIterations = 5, bool *pbStopFlag = NULL, const unsigned long nLoopKF = 0,
			const bool bRobust = true);
		void static GlobalBundleAdjustemnt(SLAM* SLAM, Map* pMap, int nIterations = 5, bool *pbStopFlag = NULL,
			const unsigned long nLoopKF = 0, const bool bRobust = true);
	};
}

#endif