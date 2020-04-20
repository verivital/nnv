load SegNet.mat;
nnvSegNet = SEGNET.parse(net, 'SetNet');
Destination = net.Connections.Destination;

