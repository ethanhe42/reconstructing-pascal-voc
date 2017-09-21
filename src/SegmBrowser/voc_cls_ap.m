function [ap, fp, tp] = voc_cls_ap(cls, scores, testset)
  VOCinit();

  VOCopts.testset = testset;
  cp=sprintf(VOCopts.annocachepath,VOCopts.testset);
  if exist(cp,'file')
    load(cp,'gtids','recs');
  else
    tic % jl added
    [gtids,t]=textread(sprintf(VOCopts.seg.imgsetpath ,VOCopts.testset),'%s %d');
    for i=1:length(gtids)
      % display progress
      if toc>1 && 0
        fprintf('%s: pr: load: %d/%d\n',cls,i,length(gtids));
        drawnow;
        tic;
      end
      
      % read annotation
      recs(i)=PASreadrecord(sprintf(VOCopts.annopath,gtids{i}));
    end
    try
      save(cp,'gtids','recs', '-V6');
    catch
      mkdir(VOCopts.annocachepath);
      save(cp,'gtids','recs', '-V6');
    end
  end

  % compile gt labels    
  gt = false(numel(gtids),1);
  cls_id = VOC09_classname_to_id(cls);
  for i=1:numel(recs)
    labels = cellfun(@VOC09_classname_to_id, {recs(i).objects.class});
    diff = [recs(i).objects.difficult];
    
    gt(i) = any(labels==cls_id & ~diff);
  end
  
  % hash image ids
  hash=VOChash_init(gtids);
  not_empty = ~cellfun(@isempty, scores);
  scores = scores(not_empty);
  ids = gtids(not_empty);
  confidence = cell2mat(scores');
   
  % map results to ground truth images
  out=ones(size(gtids))*-inf;
  tic;
  for i=1:length(ids)
      % display progress
      if toc>1 && 0
          fprintf('%s: pr: %d/%d\n',cls,i,length(ids));
          drawnow;
          tic;
      end

      % find ground truth image
      j=VOChash_lookup(hash,ids{i});
      if isempty(j)
          error('unrecognized image "%s"',ids{i});
      elseif length(j)>1
          error('multiple image "%s"',ids{i});
      else
          out(j)=confidence(i);
      end
  end

  % compute precision/recall

  [so,si]=sort(-out);
  tp=gt(si)==true;
  fp=gt(si)==false;

  cs_fp=cumsum(fp);
  cs_tp=cumsum(tp);
  rec=cs_tp/sum(gt);
  prec=cs_tp./(cs_fp+cs_tp);

  ap=VOCap(rec,prec);
  
  [v, si_back] = sort(si, 'ascend');
  tp = tp(si_back);
  fp = fp(si_back);
end