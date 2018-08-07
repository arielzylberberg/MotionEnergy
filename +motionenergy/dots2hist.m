function dotsH = dots2hist(dots,varargin)
%dotsH = dots2hist(dots,varargin)
%output: tiempo,irrelevante,relevante

N = 203;
for i=1:length(varargin)
    if isequal(varargin{i},'msize')
        N = varargin{i+1};
    end
end
nframes = length(dots); %??

% dotsHist = cell(1,nframes);
st = nan(N,N,nframes);
for k=1:nframes
    mat       = zeros(N);
    idx       = sub2ind(size(mat),ceil(dots{k}(1,:)+(N+1)/2),ceil(dots{k}(2,:)+(N+1)/2));
    mat(idx)  = 1;

%     dotsHist{k} = sparse(mat);
%     st(:,:,k) = full(dotsHist{k});
    st(:,:,k) = mat;

end

dotsH = permute(st,[3,2,1]);