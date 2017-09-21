function classname = VOC09_id_to_classname(id)
    classnames = {'aeroplane', 'bicycle', 'bird', 'boat', 'bottle', 'bus', 'car', 'cat', 'chair', ...
                             'cow', 'diningtable', 'dog', 'horse', 'motorbike', 'person', ...
                             'pottedplant', 'sheep', 'sofa', 'train', 'tvmonitor'};
       
    classname = '';
    for i=1:numel(classnames)
        if(id==i)
            classname = classnames{i};
            break;
        end
    end    
end