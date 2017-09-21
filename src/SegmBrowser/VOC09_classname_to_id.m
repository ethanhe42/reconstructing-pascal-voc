function id = VOC09_classname_to_id(classname)
    classnames = {'aeroplane', 'bicycle', 'bird', 'boat', 'bottle', 'bus', 'car', 'cat', 'chair', ...
                             'cow', 'diningtable', 'dog', 'horse', 'motorbike', 'person', ...
                             'pottedplant', 'sheep', 'sofa', 'train', 'tvmonitor'};
                         
    id = [];
    for i=1:numel(classnames)
        if(strcmp(classname , classnames{i}))
            id = i;
            break;
        end
    end
    
    if(isempty(id))
        error('no such classname');
    end
end