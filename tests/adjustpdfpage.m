%ADJUSTPDFPAGE   Adjust the figure's page size for PDF output.
%   PRINT -DPDF produces a PDF page size of A4 or letter by default.
%   Any other 'device' (i.e. -DTIFF, -DEPS, etc.) sets the bounding
%   box to match that of the figure window (or at least when
%   'PaperPositionMode' is set to 'auto'). The same is true when
%   using the figure's "Save As..." menu option.
%
%   ADJUSTPDFPAGE will set the current figure's 'PaperSize' and
%   'PaperPosition' properties such that saving the figure to PDF
%   produces the same bounding box as with other 'devices'.
%
%   ADJUSTPDFPAGE(H) will affect the figure window with handle H.
%
%   ADJUSTPDFPAGE(H,WIDTH) will set the output figure's width to
%   WIDTH centimeters. This is useful to produce graphics with a
%   consistent font size and line width.
%
%   Example:
%      h = plot(randn(100,1));
%      set(h,'linewidth',2);
%      set(gca,'fontsize',12)
%      box off
%      xlabel('repetition','fontsize',14)
%      ylabel('random value','fontsize',14)
%      adjustpdfpage(gcf,12)
%      print -dpdf figure1.pdf
%
%   Note:
%      plot(randn(100,1));
%      print -depsc2 -loose figure2.eps
%      !epstopdf figure2.eps
%   produces the same PDF file as
%      adjustpdfpage
%      print -dpdf figure2.pdf
%   (assuming you do this on a Unix computer and have epstopdf
%   installed).
%(c)2012, Cris Luengo, Centre for Image Analysis, Uppsala, Sweden.
function adjustpdfpage(figh,width)
if nargin<1
   figh = gcf;
else
   if ~isscalar(figh) || ~ishandle(figh) || ~strcmp(get(figh,'type'),'figure')
      error('First input argument should be a valid figure handle.')
   end
end
if nargin<2
   width = [];
else
   if ~isscalar(width) || ~isnumeric(width) || width<=0
      error('The WIDTH argument should be a positive scalar.')
   end
   width = double(width);
end
u = get(figh,'units');
set(figh,'units','centimeters');
p = get(figh,'position');
set(figh,'units',u);
p = p(3:4);
if ~isempty(width)
   p(2) = p(2)/p(1)*width;
   p(1) = width;
end
set(figh,'paperunits','centimeters');
set(figh,'papersize',p);
set(figh,'paperposition',[0,0,p]);