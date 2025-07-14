import torch
import numpy as np
import math
import cv2
import os
from transform import Augmentation, BaseTransform, inverse_normalize, normalize


def gen_vnnlib(im, pred, im_name):
    save_path = 'vnnlib/'
    os.makedirs(save_path, exist_ok=True)

    VOC_CLASSES = (  # always index 0
        'aeroplane', 'bicycle', 'bird', 'boat',
        'bottle', 'bus', 'car', 'cat', 'chair',
        'cow', 'diningtable', 'dog', 'horse',
        'motorbike', 'person', 'pottedplant',
        'sheep', 'sofa', 'train', 'tvmonitor')

    nn_out = pred.cpu().numpy()
    # print('NN output shape: ', nn_out.shape)
    # pred = torch.from_numpy(nn_out)
    B, abC, H, W = pred.size()
    nn_out_flat = torch.flatten(pred)

    KA = 5
    NC = 20
    prediction = pred.permute(0, 2, 3, 1).contiguous().view(B, -1, abC)
    # print(prediction.size()) # 169 x 125

    conf_pred = prediction[..., :KA].contiguous().view(B, -1, 1)
    cls_pred = prediction[..., 1 * KA: (1 + NC) * KA].contiguous().view(B, -1, NC)
    txtytwth_pred = prediction[..., (1 + NC) * KA:].contiguous().view(B, -1, 4)

    conf_pred = conf_pred[0]  # [H*W*KA, 1]
    cls_pred = cls_pred[0]  # [H*W*KA, NC]
    txtytwth_pred = txtytwth_pred[0]  # [H*W*KA, 4]

    max_conf, max_conf_idx = torch.max(conf_pred, dim=0)
    max_conf_cls = cls_pred[max_conf_idx]
    max_cls, cls_result = torch.max(max_conf_cls, dim=1)

    # print('Max_Confidence: {}'.format(max_conf.cpu().numpy()))
    # print('Max_Confidence index: {}'.format(max_conf_idx.cpu().numpy()))
    # print('Classification result: {}, {}'.format(cls_result.cpu().numpy(), VOC_CLASSES[cls_result]))

    # calculate original max confidence index
    d2 = math.floor(max_conf_idx / (W * KA))
    d3 = math.floor((max_conf_idx - d2 * W * KA) / KA)
    d1 = max_conf_idx.item() - d2 * W * KA - d3 * KA
    # print('conf idx ', d1, d2, d3)
    assert max_conf.item() == nn_out[0, d1, d2, d3]

    # calculate original max conf class index
    d4 = KA + d1 * NC + cls_result.item()
    # print('class idx ', d4, d2, d3)
    assert max_cls.item() == nn_out[0, d4, d2, d3]

    # load image
    # im = cv2.imread('dog52.jpg')  # BGR
    # im = im.transpose((2, 0, 1))[::-1]  # HWC to CHW, BGR to RGB
    # im = np.ascontiguousarray(im)  # contiguous
    # im = torch.from_numpy(im)
    iB, iC, iH, iW = im.size()
    # print('NN input shape: ', iC, iH, iW)
    # im = im.float()  # uint8 to float
    # im /= 255  # 0 - 255 to 0.0 - 1.0
    imf = torch.flatten(im)

    grid_size = H
    nclasses = NC
    boxes = KA
    epsilon = 1.0 / 255
    output_channels = boxes * (nclasses + 5)

    f = open(os.path.join(save_path, 'TinyYOLO_prop_{}_eps_1_255.vnnlib'.format(im_name)), "w")

    # Declear inputs
    # f.write("; input " + str(iC) + " x " + str(iH) + " x " + str(iW) + " brightness epsilon " + str(epsilon) + "\n\n")
    for i in range(imf.size(dim=0)):
        f.write("(declare-const X_" + str(i) + " Real)\n")

    # Declear outputs
    f.write("\n; output " + str(output_channels) + " x " + str(grid_size) + " x " + str(grid_size) + "\n\n")
    for i in range(nn_out_flat.size(dim=0)):
        f.write("(declare-const Y_" + str(i) + " Real)\n")

    # Declear input constraints
    f.write("\n; input constraints\n\n")
    perturb = torch.ones_like(im) * epsilon
    # print(perturb)
    x_nat = inverse_normalize(im.detach().clone())
    im_ub = normalize(torch.max(torch.min(x_nat + perturb, torch.ones_like(im)), torch.zeros_like(im)))
    im_lb = normalize(torch.max(torch.min(x_nat - perturb, torch.ones_like(im)), torch.zeros_like(im)))
    # print(x_net + perturb)
    imf_ub = torch.flatten(im_ub)
    imf_lb = torch.flatten(im_lb)
    assert torch.all(imf_ub >= imf_lb)

    for i in range(imf.size(dim=0)):
        # val = imf[i].item()
        # ub = min(val + epsilon, 1.0)
        # lb = max(val - epsilon, 0.0)
        ub = imf_ub[i].item()
        lb = imf_lb[i].item()
        f.write("(assert (<= X_" + str(i) + " " + str(ub) + "))\n")
        f.write("(assert (>= X_" + str(i) + " " + str(lb) + "))\n")

    # Declear output constraints: only check one grid cell, either no object detected or classification changes
    f.write("\n; output constraints\n\n")

    x = d3
    y = d2
    # Property: the grid cell containnig the center of the object will not change
    # f.write("(assert (and\n")
    # confidence constraints: one of the bounding boxes has a confiden score > threshold
    threshold = -1 # signmoid(threshold) = 0.1
    for idx in range(boxes):
      # f.write("  (or\n")
      conf_idx = idx*grid_size*grid_size + y*grid_size + x
      f.write("(assert (<= Y_" + str(conf_idx) + " " + str(threshold) + "))\n")
    #   # classification constraints: the classification does not change
    #   max_class_idx = (boxes+idx*nclasses+cls_result.item())*grid_size*grid_size + y*grid_size + x
    #   for ic in range(boxes+idx*nclasses, boxes+(idx+1)*nclasses):
    #      class_idx = ic*grid_size*grid_size + y*grid_size + x
    #      if class_idx != max_class_idx:
    #         f.write("     (and (<= Y_" + str(max_class_idx) + " " + "Y_" + str(class_idx) + "))\n")
    #   f.write("  )\n")
    # f.write("))\n")
    f.close()