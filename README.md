
# Cy-RoPS

This code implements localization algorithm in underground garages with Cy-RoPS feature.

For more detailed information please refer to the following paper:

	Zhongxing Tao, Jianru Xue, Di Wang, Shuyang Zhang, Dixiao Cui, and Shaoyi Du, "Accurate Localization in Underground Garages via Cylinder Feature based Map Matching,” in Proc. IEEE Intelligent Vehicles Symposium (IV), 2018.

In this paper, we propose a cylinder rotational projection statistics feature (Cy-RoPS), which is motivated by local surface feature based 3D object recognition method - RoPS presented by Yulan Guo:

	Y. Guo, F. A. Sohel, M. Bennamoun, M. Lu, and J. Wan, “Rotational projection statistics for 3d local surface description and object recognition,”International Journal of Computer Vision, vol. 105, no. 1, pp. 63–86, 2013.

# Function in Matlab:

	MainTest.m		-	Localizaton demo

	GrdSegmentFun.m		-	Ground segmentation of Lidar data

	RecoCylinderFun.m	-	Cy-RoPS function
	
	IcpFun.m		-	Improved point to plane ICP algorithm, no need to create a KD tree every frame

	StatisticalRatioFun.m	- 	Calculate the ratio which approximates the rate of matching

	RecoAndDetectDemo.m	-	Demo for detecting pillars and wall in an underground garage.

#  Affiliation

Emails : newwaves@163.com

Laboratory of Visual Cognitive Computing and Intelligent Vehicle

Institute of Artificial Intelligence and Robotics

Xi'an Jiaotong University

710049 Xi'an, Shaanxi, 710049, P.R. China
