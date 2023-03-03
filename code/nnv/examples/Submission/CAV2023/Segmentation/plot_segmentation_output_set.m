function plot_segmentation_output_set(obj, ind, save_name)   
    
    % Plot verification result of segmentation task
    if isempty(obj.verifiedOutputSet)
        error("Verified Output Reachable set is empty, please perform verify method first");
    end
    
    ps = obj.pixelClassificationReachSet;
    rs = obj.verifiedOutputSet;
    gr = obj.groundTruthSegIms;
    
    if ind > length(rs) || ind < 1
        error("Invalid index");
    end
    
    pl_rs = rs{ind};
    pl_gr = gr{ind};         
    px_rs = ps{ind};
    gr_RGB = label2rgb(pl_gr);
    rs_RGB = label2rgb(pl_rs);
    px_RGB = label2rgb(px_rs);
    
    gr_unique = unique(pl_gr);
    m1 = length(gr_unique);
    [IND1,in_map1] = rgb2ind(gr_RGB, m1);
    map1 = obj.getColorMap(pl_gr, IND1, in_map1);
    classes1 = obj.getClasses(gr_unique); 
    
    rs_unique = unique(pl_rs);
    m2 = length(rs_unique);
    [IND2,in_map2] = rgb2ind(rs_RGB, m2);
    map2 = obj.getColorMap(pl_rs, IND2, in_map2);
    classes2 = obj.getClasses(rs_unique); 
    
    px_unique = unique(px_rs);
    m3 = length(px_unique);
    [IND3,in_map3] = rgb2ind(px_RGB, m3);
    map3 = obj.getColorMap(px_rs, IND3, in_map3);
    classes3 = obj.getClasses(px_unique); 

    % Set figure size before saving
    figSize = [400 400];
    
    f1 = figure;
%     ax1 = subplot(1,3,1);
    ax1 = gca;
    imshow(gr_RGB);
    colormap(ax1,map1);
    cbh1 = colorbar(ax1);
    xtick = 1/(2*m1):1/m1:1;
    cbh1.Ticks = xtick;               
    cbh1.TickLabels = classes1;
%     title("Input image without attack");
    truesize(figSize);
    exportgraphics(f1, save_name+"_1.pdf", 'ContentType', 'vector');

%     ax3 = subplot(1,3,2);
    f2 = figure;
    ax3 = gca;
    imshow(px_RGB);
    colormap(ax3,map3);
    cbh1 = colorbar(ax3);
    xtick = 1/(2*m3):1/m3:1;
    cbh1.Ticks = xtick;               
    cbh1.TickLabels = classes3;
%     title("Computed Reach Set");
    truesize(figSize);
    exportgraphics(f2, save_name+"_2.pdf", 'ContentType', 'vector');
    
%     ax2 = subplot(1,3,3);
    f3 = figure;
    ax2 = gca;
    imshow(rs_RGB);
    colormap(ax2,map2);
    cbh2 = colorbar(ax2);
    xtick = 1/(2*m2):1/m2:1;
    cbh2.Ticks = xtick;               
    cbh2.TickLabels = classes2;
%     title("Verified reach set");
    truesize(figSize);
    exportgraphics(f3, save_name+"_3.pdf", 'ContentType', 'vector');
    
end

