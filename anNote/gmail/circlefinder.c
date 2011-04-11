#include "cv.h"
#include "highgui.h"

#include <stdio.h>

//
// We need this to be high enough to get rid of things that are too small too
// have a definite shape.  Otherwise, they will end up as ellipse false positives.
//
#define MIN_AREA 100.00    
//
// One way to tell if an object is an ellipse is to look at the relationship
// of its area to its dimensions.  If its actual occupied area can be estimated
// using the well-known area formula Area = PI*A*B, then it has a good chance of
// being an ellipse.
//
// This value is the maximum permissible error between actual and estimated area.
//
#define MAX_TOL  100.00

int main( int argc, char** argv )
{
    IplImage* src;
    // the first command line parameter must be file name of binary (black-n-white) image
    if( argc == 2 && (src=cvLoadImage(argv[1], 0))!= 0)
    {
        IplImage* dst  = cvCreateImage( cvGetSize(src), 8, 3 );
        CvMemStorage* storage = cvCreateMemStorage(0);
        CvSeq* contour = 0;    
        cvThreshold( src, src, 1, 255, CV_THRESH_BINARY );
        //
        // Invert the image such that white is foreground, black is background.
        // Dilate to get rid of noise.
        //
        cvXorS(src, cvScalar(255), src);
        cvDilate(src, src, NULL, 2);    
        cvFindContours( src, storage, &contour, sizeof(CvContour), CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE );
        cvZero( dst );

        for( ; contour != 0; contour = contour->h_next )
        {
            double actual_area = fabs(cvContourArea(contour));
            if (actual_area < MIN_AREA)
                continue;

            //
            // FIXME:
            // Assuming the axes of the ellipse are vertical/perpendicular.
            //
            CvRect rect = ((CvContour *)contour)->rect;
            int A = rect.width / 2; 
            int B = rect.height / 2;
            double estimated_area = M_PI * A * B;
            double error = fabs(actual_area - estimated_area);    
            if (error > MAX_TOL)
                continue;    
            printf
            (
                 "center x: %d y: %d A: %d B: %d\n",
                 rect.x + A,
                 rect.y + B,
                 A,
                 B
            );

            CvScalar color = CV_RGB( rand() % 255, rand() % 255, rand() % 255 );
            cvDrawContours( dst, contour, color, color, -1, CV_FILLED, 8 );
        }

        cvSaveImage("coins.png", dst);
    }
}
