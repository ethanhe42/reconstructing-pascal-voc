function kp_struct = read_brkl_kp(filename)
  theStruct = parseXML(filename);
  
  for i=1:numel(theStruct.Children)  
      if(strcmp(theStruct.Children(i).Name, 'image'))
          kp_struct.img_name = theStruct.Children(i).Children(1).Data;
      elseif(strcmp(theStruct.Children(i).Name, 'voc_id'))
          kp_struct.voc_id = str2double(theStruct.Children(i).Children.Data);
      elseif(strcmp(theStruct.Children(i).Name, 'visible_bounds'))
            vb = theStruct.Children(i);
            if(~isempty(vb.Attributes))
                kp_struct.visible_bounds = [str2double(vb.Attributes(3).Value) ...
                    str2double(vb.Attributes(4).Value) ...
                    str2double(vb.Attributes(1).Value) ...
                    str2double(vb.Attributes(2).Value)];
                
                assert(strcmp(vb.Attributes(3).Name, 'xmin'));
                assert(strcmp(vb.Attributes(4).Name, 'ymin'));
            else
                kp_struct.visible_bounds = [];
            end
      elseif(strcmp(theStruct.Children(i).Name, 'subcategory'))
          kp_struct.subclass = theStruct.Children(i).Children.Data;          
          id_kp = 12;
      elseif(strcmp(theStruct.Children(i).Name, 'category'))
          kp_struct.class = theStruct.Children(i).Children.Data;
      elseif(strcmp(theStruct.Children(i).Name, 'keypoints'))
          kp = theStruct.Children(i);
          to_remove = false(numel(kp.Children),1);
          for j=1:numel(kp.Children)
              if(strcmp(kp.Children(j).Name, '#text'))
                  to_remove(j) = true;
              end
          end
          
          kp.Children(to_remove) = [];
          kp_struct.kp = zeros(3,numel(kp.Children));
          kp_struct.kp_name = cell(numel(kp.Children),1);
          for j=1:numel(kp.Children)
              kp_struct.kp_coords(:,j) = [str2double(kp.Children(j).Attributes(3).Value) ...
                  str2double(kp.Children(j).Attributes(4).Value) ...
                  str2double(kp.Children(j).Attributes(5).Value)]';
              kp_struct.kp_name{j} = kp.Children(j).Attributes(1).Value;
              kp_struct.visib(j) = str2double(kp.Children(j).Attributes(2).Value);
          end
      end
  end    
end