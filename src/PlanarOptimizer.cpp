#include <PlanarOptimizer.h>
#include <SLAM.h>
#include <Frame.h>
#include <KeyFrame.h>
#include <Map.h>
#include <MapPoint.h>
#include <Converter.h>
#include <LoopCloser.h>
#include <LabelInfo.h>

#include "g2o/core/block_solver.h"
#include <g2o/solvers/csparse/linear_solver_csparse.h>
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/types/sba/types_six_dof_expmap.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/types/sim3/types_seven_dof_expmap.h"

#include <PlaneNode.h>
#include "PlaneEstimator.h"
#include "SemanticProcessor.h"

namespace EdgeSLAM {
	void PlanarOptimizer::GlobalBundleAdjustemnt(SLAM* SLAM, Map* pMap, int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
	{
		std::vector<KeyFrame*> vpKFs = pMap->GetAllKeyFrames();
		std::vector<MapPoint*> vpMP = pMap->GetAllMapPoints();
		std::cout << "GBA::Start" << std::endl;
		BundleAdjustment(SLAM, vpKFs, vpMP, nIterations, pbStopFlag, nLoopKF, bRobust);
		std::cout << "GBA::End" << std::endl;
	}


	void PlanarOptimizer::BundleAdjustment(SLAM* SLAM, const std::vector<KeyFrame *> &vpKFs, const std::vector<MapPoint *> &vpMP,
		int nIterations, bool* pbStopFlag, const unsigned long nLoopKF, const bool bRobust)
	{
		std::vector<bool> vbNotIncludedMP;
		vbNotIncludedMP.resize(vpMP.size());

		g2o::SparseOptimizer optimizer;
		auto linear_solver = std::make_unique<g2o::LinearSolverEigen<g2o::BlockSolver_6_3::PoseMatrixType>>();
		auto block_solver = std::make_unique<g2o::BlockSolver_6_3>(std::move(linear_solver));
		auto algorithm = new g2o::OptimizationAlgorithmLevenberg(std::move(block_solver));
		optimizer.setAlgorithm(algorithm);

		if (pbStopFlag)
			optimizer.setForceStopFlag(pbStopFlag);

		long unsigned int maxKFid = 0;
		

		// Set KeyFrame vertices
		for (size_t i = 0; i<vpKFs.size(); i++)
		{
			KeyFrame* pKF = vpKFs[i];
			if (pKF->isBad())
				continue;
			g2o::VertexSE3Expmap * vSE3 = new g2o::VertexSE3Expmap();
			vSE3->setEstimate(Converter::toSE3Quat(pKF->GetPose()));
			vSE3->setId(pKF->mnId);
			vSE3->setFixed(pKF->mnId == 0);
			optimizer.addVertex(vSE3);
			if (pKF->mnId>maxKFid)
				maxKFid = pKF->mnId;
		}

		//Set Plane vertices
		long unsigned int maxPlaneid = 0;
		long unsigned int floorID = 0;
		long unsigned int ceilID = 0;

		std::list<SemanticSLAM::Plane*> lLocalPlanes;
		std::list<SemanticSLAM::PlaneEstRes*> lLocalWallEst;
		std::map<SemanticSLAM::PlaneEstRes*, SemanticSLAM::Plane*> mapWallPlanes;
		if (SemanticSLAM::PlaneEstimator::GlobalFloor->nScore > 0) {
			lLocalPlanes.push_back(SemanticSLAM::PlaneEstimator::GlobalFloor);
			floorID = SemanticSLAM::PlaneEstimator::GlobalFloor->mnID;

			////벽
			auto mapWallDatas = SemanticSLAM::PlaneEstimator::GlobalNormalMPs.Get();
			std::vector<std::pair<int, SemanticSLAM::PlaneEstRes*>> vec(mapWallDatas.begin(), mapWallDatas.end());
			std::sort(vec.begin(), vec.end(), SemanticSLAM::PlaneEstimator::cmp);
			cv::Mat Rsp = SemanticSLAM::PlaneEstimator::CalcPlaneRotationMatrix(SemanticSLAM::PlaneEstimator::GlobalFloor->param).clone();

			float max_dist = -FLT_MAX;
			for (int i = 0; i < vec.size(); i++) {
				auto p = vec[i].second;
				bool bwall = true;

				auto tempWall = p;
				if(tempWall->nData > 0)
					std::cout << "wall data = " << tempWall->nData << std::endl;
				while(bwall)
				{
					if (tempWall->nData < 300) {
						bwall = false;
						break;
					}

					tempWall->SetData();
					tempWall->ConvertWallData(Rsp);
					SemanticSLAM::PlaneEstRes* tempIn = new SemanticSLAM::PlaneEstRes();
					SemanticSLAM::PlaneEstRes* tempOut = new SemanticSLAM::PlaneEstRes();

					bool bres = SemanticSLAM::PlaneEstimator::PlaneInitialization(tempWall, tempIn, tempOut, 1500, 0.01, 0.2);
					if (bres) {
						cv::Mat param = cv::Mat::zeros(4, 1, CV_32FC1);
						param.at<float>(0) = tempIn->param.at<float>(0);
						param.at<float>(2) = tempIn->param.at<float>(1);
						param.at<float>(3) = tempIn->param.at<float>(2);

						SemanticSLAM::Plane* pa = new SemanticSLAM::Plane();
						pa->type = SemanticSLAM::PlaneType::WALL;
						pa->param = param.clone();
						pa->normal = pa->param.rowRange(0, 3);
						pa->dist = cv::norm(pa->normal);
						pa->nScore = tempIn->data.rows;
						pa->count = 1;
						tempIn->nPlaneID = pa->mnID;

						//최적화시 정보 추가
						lLocalWallEst.push_back(tempIn);
						lLocalPlanes.push_back(pa);
						mapWallPlanes[tempIn] = pa;
						tempWall = tempOut;

						float dist = param.at<float>(3);
						if (dist > max_dist) {
							max_dist = dist;
							SemanticSLAM::PlaneEstimator::MainWall = pa;
						}
						/*if (!SemanticSLAM::PlaneEstimator::MainWall) {
							
						}*/
					}
					else {
						bwall = false;
						break;
					}
				}
			}//for
		}
		if (SemanticSLAM::PlaneEstimator::GlobalCeil->nScore > 0) {
			lLocalPlanes.push_back(SemanticSLAM::PlaneEstimator::GlobalCeil);
			ceilID = SemanticSLAM::PlaneEstimator::GlobalCeil->mnID;
		}
		for (auto lit = lLocalPlanes.begin(), lend = lLocalPlanes.end(); lit != lend; lit++) {
			SemanticSLAM::Plane* pPlane = *lit;
			g2o::PlaneVertex* vPlane = new g2o::PlaneVertex();
			vPlane->setEstimate(Converter::toVector6d(pPlane->param));
			vPlane->setId(maxKFid + pPlane->mnID);
			vPlane->setFixed(true);
			optimizer.addVertex(vPlane);
			if (pPlane->mnID > maxPlaneid)
				maxPlaneid = pPlane->mnID;
		}

		const float thHuber2D = sqrt(5.99);

		// Set MapPoint vertices
		std::vector<MapPoint*> vpPlaneEdgeMP;
		std::vector<int> vpPlaneEdgdID;
		std::vector<g2o::PlaneBAEdge*> vpPlaneEdge;
		std::vector<SemanticSLAM::Plane*> vpPlanes;
		for (size_t i = 0; i<vpMP.size(); i++)
		{
			MapPoint* pMP = vpMP[i];
			if (pMP->isBad())
				continue;
			g2o::VertexPointXYZ* vPoint = new g2o::VertexPointXYZ();
			vPoint->setEstimate(Converter::toVector3d(pMP->GetWorldPos()));
			const int id = pMP->mnId + maxKFid + maxPlaneid + 1;
			vPoint->setId(id);
			vPoint->setMarginalized(true);
			optimizer.addVertex(vPoint);

			const std::map<KeyFrame*, size_t> observations = pMP->GetObservations();

			int nEdges = 0;
			//SET EDGES
			for (std::map<KeyFrame*, size_t>::const_iterator mit = observations.begin(); mit != observations.end(); mit++)
			{

				KeyFrame* pKF = mit->first;
				if (pKF->isBad() || pKF->mnId>maxKFid)
					continue;

				nEdges++;

				const cv::KeyPoint &kpUn = pKF->mvKeysUn[mit->second];

				Eigen::Matrix<double, 2, 1> obs;
				obs << kpUn.pt.x, kpUn.pt.y;

				g2o::EdgeSE3ProjectXYZ* e = new g2o::EdgeSE3ProjectXYZ();

				e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
				e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pKF->mnId)));
				e->setMeasurement(obs);
				const float &invSigma2 = pKF->mvInvLevelSigma2[kpUn.octave];
				e->setInformation(Eigen::Matrix2d::Identity()*invSigma2);

				if (bRobust)
				{
					g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
					e->setRobustKernel(rk);
					rk->setDelta(thHuber2D);
				}

				e->fx = pKF->fx;
				e->fy = pKF->fy;
				e->cx = pKF->cx;
				e->cy = pKF->cy;

				optimizer.addEdge(e);
			}

			if (nEdges == 0)
			{
				optimizer.removeVertex(vPoint);
				vbNotIncludedMP[i] = true;
			}
			else
			{
				vbNotIncludedMP[i] = false;
			}

			//if (floorID>0 && pMP->mnLabelID == (int)SemanticSLAM::PlaneType::FLOOR) {
			//	if (vbNotIncludedMP[i])
			//		continue;
			//	g2o::PlaneBAEdge* e = new g2o::PlaneBAEdge();
			//	e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
			//	e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(floorID + maxKFid)));
			//	e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());

			//	g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
			//	//e->setRobustKernel(rk);
			//	//rk->setDelta(thPHuberPlane);

			//	optimizer.addEdge(e);
			//	vpPlaneEdge.push_back(e);
			//	vpPlaneEdgeMP.push_back(pMP);
			//	vpPlaneEdgePlane.push_back(SemanticSLAM::PlaneEstimator::GlobalFloor);
			//}
		}

		/////floor mp
		{
			auto vecFloorMPs = SemanticSLAM::PlaneEstimator::FloorAllData->vecMPs;
			for (size_t i = 0; i < vecFloorMPs.size(); i++)
			{
				MapPoint* pMP = vecFloorMPs[i];
				if (pMP->isBad())
					continue;
				const int id = pMP->mnId + maxKFid + maxPlaneid + 1;
				g2o::PlaneBAEdge* e = new g2o::PlaneBAEdge();
				e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
				e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(floorID + maxKFid)));
				e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());

				g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
				//e->setRobustKernel(rk);
				//rk->setDelta(thPHuberPlane);

				optimizer.addEdge(e);
				vpPlaneEdge.push_back(e);
				vpPlaneEdgeMP.push_back(pMP);
				vpPlanes.push_back(SemanticSLAM::PlaneEstimator::GlobalFloor);
				vpPlaneEdgdID.push_back(SemanticSLAM::PlaneEstimator::GlobalFloor->mnID);
			}
		}

		{
			auto vecCeilMPs = SemanticSLAM::PlaneEstimator::CeilAllData->vecMPs;
			for (size_t i = 0; i < vecCeilMPs.size(); i++)
			{
				MapPoint* pMP = vecCeilMPs[i];
				if (pMP->isBad())
					continue;
				const int id = pMP->mnId + maxKFid + maxPlaneid + 1;
				g2o::PlaneBAEdge* e = new g2o::PlaneBAEdge();
				e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
				e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(ceilID + maxKFid)));
				e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());

				g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
				//e->setRobustKernel(rk);
				//rk->setDelta(thPHuberPlane);

				optimizer.addEdge(e);
				vpPlaneEdge.push_back(e);
				vpPlaneEdgeMP.push_back(pMP);
				vpPlanes.push_back(SemanticSLAM::PlaneEstimator::GlobalCeil);
				vpPlaneEdgdID.push_back(SemanticSLAM::PlaneEstimator::GlobalCeil->mnID);
			}

		}

		{
			for (auto lter = lLocalWallEst.begin(), lend = lLocalWallEst.end(); lter != lend; lter++) {
				auto pres = *lter;
				auto vecMPs = pres->vecMPs;
				auto plane = mapWallPlanes[pres];

				for (size_t i = 0; i < vecMPs.size(); i++)
				{
					MapPoint* pMP = vecMPs[i];
					if (pMP->isBad())
						continue;
					const int id = pMP->mnId + maxKFid + maxPlaneid + 1;
					g2o::PlaneBAEdge* e = new g2o::PlaneBAEdge();
					e->setVertex(1, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(id)));
					e->setVertex(0, dynamic_cast<g2o::OptimizableGraph::Vertex*>(optimizer.vertex(pres->nPlaneID + maxKFid)));
					e->setInformation(Eigen::Matrix<double, 1, 1>::Identity());

					g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
					//e->setRobustKernel(rk);
					//rk->setDelta(thPHuberPlane);

					optimizer.addEdge(e);
					vpPlaneEdge.push_back(e);
					vpPlaneEdgeMP.push_back(pMP);
					vpPlaneEdgdID.push_back(pres->nPlaneID);
					vpPlanes.push_back(plane);
				}
			}
		}


		// Optimize!
		optimizer.initializeOptimization();
		optimizer.optimize(nIterations);

		// Recover optimized data

		//Keyframes
		for (size_t i = 0; i<vpKFs.size(); i++)
		{
			KeyFrame* pKF = vpKFs[i];
			if (pKF->isBad())
				continue;
			g2o::VertexSE3Expmap* vSE3 = static_cast<g2o::VertexSE3Expmap*>(optimizer.vertex(pKF->mnId));
			g2o::SE3Quat SE3quat = vSE3->estimate();
			if (nLoopKF == 0)
			{
				pKF->SetPose(Converter::toCvMat(SE3quat));
			}
			else
			{
				pKF->mTcwGBA.create(4, 4, CV_32F);
				Converter::toCvMat(SE3quat).copyTo(pKF->mTcwGBA);
				pKF->mnBAGlobalForKF = nLoopKF;
			}
		}

		/*for (auto iter = lLocalPlanes.begin(), iend = lLocalPlanes.end(); iter != iend; iter++) {
			auto plane = *iter;
			g2o::PlaneVertex* vPlane = static_cast<g2o::PlaneVertex*>(optimizer.vertex(maxKFid + plane->mnID));
			auto param = Converter::toCvMat(vPlane->estimate());
			if (plane->type == SemanticSLAM::PlaneType::FLOOR) {
				std::cout << "FLOOR==" << plane->param.t() << " " << param.t() << std::endl;
			}
			if (plane->type == SemanticSLAM::PlaneType::CEIL) {
				std::cout << "CEIL==" << plane->param.t() << " " << param.t() << std::endl;
			}
			plane->param = param.rowRange(0, 4);
			plane->normal = plane->param.rowRange(0, 3);
			plane->dist = cv::norm(plane->normal);
		}*/
		
		//floor points
		int nFloorError = 0;
		/*if (SemanticSLAM::PlaneEstimator::GlobalFloor) {
		g2o::PlaneVertex* vPlane = static_cast<g2o::PlaneVertex*>(optimizer.vertex( maxKFid + floorID));
		std::cout<<Converter::toCvMat(vPlane->estimate()).rowRange(0, 4).t()<< SemanticSLAM::PlaneEstimator::GlobalFloor->param.t();
		}*/

		float pthresh = 0.001;//0.0001;
		std::map<int, cv::Mat> mapDatas, mapDatas2, mapDatas3;
		std::map<SemanticSLAM::Plane*, std::vector<MapPoint*>> mapWallMPs;
		for (size_t i = 0, iend = vpPlaneEdge.size(); i < iend; i++) {
			g2o::PlaneBAEdge* e = vpPlaneEdge[i];
			MapPoint* pMP = vpPlaneEdgeMP[i];
			auto plane = vpPlanes[i];
			if (!pMP)
				continue;
			if (pMP->isBad())
				continue;
			if (e->chi2() > pthresh) {
				nFloorError++;
				//mapDatas2[pMP->mnId] = pMP->GetWorldPos();
			}
			else {
				//시각화
				if (vpPlaneEdgdID[i] == floorID) {
					pMP->mnPlaneID = (int)StructureLabel::FLOOR;
					//mapDatas[pMP->mnId] = pMP->GetWorldPos();
				}else if (vpPlaneEdgdID[i] == ceilID) {
					pMP->mnPlaneID = (int)StructureLabel::CEIL;
					//mapDatas2[pMP->mnId] = pMP->GetWorldPos();
				}
				else {
					pMP->mnPlaneID = (int)StructureLabel::WALL;
					//mapDatas3[pMP->mnId] = pMP->GetWorldPos();

					mapWallMPs[plane].push_back(pMP);
					//해당 플레인에 대한 인접한 키프레임.
					//맵포인트에 옵저베이션으로부터 시작하기.
					//이게 어떤 평면에 해당하는지, 아이디가 무엇인지 알 수 있ㅇ므.
				}
			}
		}

		std::cout << "Calculate Connected KF and Wall Planes" << std::endl;
		////KF 연결하기
		//같은 평면에 속하는 MP들
		std::map<KeyFrame*, std::set<SemanticSLAM::Plane*>> mapKeyFrameNPlanes;
		for (auto iter = mapWallMPs.begin(), iend = mapWallMPs.end(); iter != iend; iter++) {
			std::map<KeyFrame*, int> keyframeCounter;
			std::vector<KeyFrame*> vpLocalKFs;
			std::set<KeyFrame*> spLocalKFs;
			
			auto vec = iter->second;
			auto plane = iter->first;

			for (int j = 0; j < vec.size(); j++) {
				auto pMP = vec[j];
				if (!pMP)
					continue;
				if (pMP->isBad())
					continue;
				const std::map<KeyFrame*, size_t> observations = pMP->GetObservations();
				for (std::map<KeyFrame*, size_t>::const_iterator it = observations.begin(), itend = observations.end(); it != itend; it++)
					keyframeCounter[it->first]++;

				for (std::map<KeyFrame*, int>::const_iterator it = keyframeCounter.begin(), itEnd = keyframeCounter.end(); it != itEnd; it++)
				{
					KeyFrame* pKF = it->first;
					if (pKF->isBad())
						continue;
					spLocalKFs.insert(it->first);
					vpLocalKFs.push_back(it->first);

				}

				//for (size_t k = 0, kend = vpLocalKFs.size(); k < kend; k++)
				//{
				//	KeyFrame* pKF = vpLocalKFs[k];// *itKF;

				//	const std::vector<KeyFrame*> vNeighs = pKF->GetBestCovisibilityKeyFrames(10);

				//	for (std::vector<KeyFrame*>::const_iterator itNeighKF = vNeighs.begin(), itEndNeighKF = vNeighs.end(); itNeighKF != itEndNeighKF; itNeighKF++)
				//	{
				//		KeyFrame* pNeighKF = *itNeighKF;
				//		if (pNeighKF && !pNeighKF->isBad() && !spLocalKFs.count(pNeighKF))
				//		{
				//			spLocalKFs.insert(pNeighKF);
				//			continue;
				//		}
				//	}
				//}//vpLocalKFs

			}//for mp

			//현재 벽과 연결된 MP와 연결된 키프레임들에 벽을 연결
			for (auto iter = spLocalKFs.begin(), iend = spLocalKFs.end(); iter != iend; iter++) {
				auto pKF = *iter;
				if(!mapKeyFrameNPlanes[pKF].count(plane))
					mapKeyFrameNPlanes[pKF].insert(plane);
			}
		}//for wall

		std::cout << "Update Connected KF and Wall Planes = " << mapWallPlanes.size()<<" "<< mapWallMPs.size() << std::endl;
		//키프레임으로부터 평면 연결해야 함.
		for (auto iter = mapKeyFrameNPlanes.begin(), iend = mapKeyFrameNPlanes.end(); iter != iend; iter++) {
			auto pKF = iter->first;
			auto setPlanes = iter->second;
			SemanticSLAM::PlaneEstimator::mPlaneConnections.Update(pKF, setPlanes);
		}
		
		/*if (mapDatas.size() > 0) {
			SLAM->TemporalDatas2.Update("GBAFloor", mapDatas);
		}
		if (mapDatas2.size() > 0) {
			SLAM->TemporalDatas2.Update("GBACeil", mapDatas2);
		}
		if (mapDatas3.size() > 0) {
			SLAM->TemporalDatas2.Update("GBAWall", mapDatas3);
		}
		std::cout << "FloorMPs = " << vpPlaneEdgeMP.size() <<" "<<nFloorError<< std::endl;*/

		//Points
		for (size_t i = 0; i<vpMP.size(); i++)
		{
			if (vbNotIncludedMP[i])
				continue;

			MapPoint* pMP = vpMP[i];

			if (pMP->isBad())
				continue;
			g2o::VertexPointXYZ* vPoint = static_cast<g2o::VertexPointXYZ*>(optimizer.vertex(pMP->mnId + maxKFid + maxPlaneid + 1));

			if (nLoopKF == 0)
			{
				pMP->SetWorldPos(Converter::toCvMat(vPoint->estimate()));
				pMP->UpdateNormalAndDepth();
			}
			else
			{
				pMP->mPosGBA.create(3, 1, CV_32F);
				Converter::toCvMat(vPoint->estimate()).copyTo(pMP->mPosGBA);
				pMP->mnBAGlobalForKF = nLoopKF;
			}
		}

		//////평면
		//std::map<int, cv::Mat> mapDatas4;
		//auto mapWallDatas = SemanticSLAM::PlaneEstimator::GlobalNormalMPs.Get();
		//std::vector<std::pair<int, SemanticSLAM::PlaneEstRes*>> vec(mapWallDatas.begin(), mapWallDatas.end());
		//std::sort(vec.begin(), vec.end(), SemanticSLAM::PlaneEstimator::cmp);

		//cv::Mat Rsp = SemanticSLAM::PlaneEstimator::CalcPlaneRotationMatrix(SemanticSLAM::PlaneEstimator::GlobalFloor->param).clone();

		//for (int i = 0; i < vec.size(); i++) {
		//	auto p = vec[i].second;
		//	if (p->nData > 1000) {
		//		p->SetData();
		//		p->ConvertWallData(Rsp);
		//		SemanticSLAM::PlaneEstRes* tempIn = new SemanticSLAM::PlaneEstRes();
		//		SemanticSLAM::PlaneEstRes* tempOut = new SemanticSLAM::PlaneEstRes();
		//		bool bres = SemanticSLAM::PlaneEstimator::PlaneInitialization(p, tempIn, tempOut, 1500, 0.01, 0.2);
		//		if (bres) {
		//			auto vecMPs = tempIn->vecMPs;
		//			for (size_t i = 0, iend = vecMPs.size(); i < iend; i++) {
		//				auto pMP = vecMPs[i];
		//				if (pMP->isBad())
		//					continue;
		//				if (pMP->mnPlaneCount > 10)
		//					mapDatas4[pMP->mnId] = pMP->GetWorldPos();
		//			}
		//		}
		//		delete tempIn;
		//		delete tempOut;
		//		//std::cout << p->nData<<std::endl;
		//	}//if n
		//}//for
		//SLAM->TemporalDatas2.Update("GBAWall", mapDatas4);

	}//BA

}