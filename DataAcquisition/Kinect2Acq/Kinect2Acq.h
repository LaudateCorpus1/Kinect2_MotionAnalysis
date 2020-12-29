//------------------------------------------------------------------------------
// <copyright file="Kinect2Acq.h" company="Microsoft">
//     The Kinect for Windows APIs used here are preliminary and subject to change
//     Copyright (c) Microsoft Corporation.  All rights reserved.
// </copyright>
//------------------------------------------------------------------------------

#pragma once

#include "resource.h"
#include <stdio.h>
#include <stdlib.h>

class CKinect2Acq
{
    //Kinect Sensor parameters
	static const int        cDepthWidth  = 512;
    static const int        cDepthHeight = 424;

	static const int        cColorWidth  = 1920;
	static const int        cColorHeight = 1080;

	//Video file parameters
	const UINT32 VIDEO_WIDTH           = cColorWidth;
	const UINT32 VIDEO_HEIGHT          = cColorHeight;
	const UINT32 VIDEO_FPS             = 30;							//frame rate of device
	const UINT64 VIDEO_FRAME_DURATION  = 10 * 1000 * 1000 / VIDEO_FPS;	//duration of frame in 100 nanosecond increments
	const UINT32 VIDEO_BIT_RATE        = 800000;
	const GUID   VIDEO_ENCODING_FORMAT = MFVideoFormat_H264;			//video encoding format
	const GUID   VIDEO_INPUT_FORMAT    = MFVideoFormat_RGB32;
	const UINT32 VIDEO_PELS            = VIDEO_WIDTH * VIDEO_HEIGHT;


public:
    /// <summary>
    /// Constructor
    /// </summary>
    CKinect2Acq();

    /// <summary>
    /// Destructor
    /// </summary>
    ~CKinect2Acq();

    /// <summary>
    /// Handles window messages, passes most to the class instance to handle
    /// </summary>
    /// <param name="hWnd">window message is for</param>
    /// <param name="uMsg">message</param>
    /// <param name="wParam">message data</param>
    /// <param name="lParam">additional message data</param>
    /// <returns>result of message processing</returns>
    static LRESULT CALLBACK MessageRouter(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

    /// <summary>
    /// Handle windows messages for a class instance
    /// </summary>
    /// <param name="hWnd">window message is for</param>
    /// <param name="uMsg">message</param>
    /// <param name="wParam">message data</param>
    /// <param name="lParam">additional message data</param>
    /// <returns>result of message processing</returns>
    LRESULT CALLBACK        DlgProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
	
    /// <summary>
    /// Creates the main window and begins processing
    /// </summary>
    /// <param name="hInstance"></param>
    /// <param name="nCmdShow"></param>
    int                     Run(HINSTANCE hInstance, int nCmdShow);

private:

	//basic helper variables
    HWND                    m_hWnd;
    INT64                   m_nStartTime;
    INT64                   m_nLastCounter;
    static double           m_fFreq;
	//static float            m_PulseLength;
    DWORD                   m_nNextStatusTime;
    DWORD                   m_nFramesSinceUpdate;
	LONG					m_sourceStride;
	bool                    m_bSaveScreenshot;

    // Current Kinect
    IKinectSensor*          m_pKinectSensor;
    ICoordinateMapper*      m_pCoordinateMapper;

    // Frame reader
	IMultiSourceFrameReader*m_pMultiSourceFrameReader;

	//File Writing
	bool                    m_bRecordData;
	FILE*					m_pDataFile;		//file for writing out camera xyz coords to binary file
	FILE*					m_pImageCoordsFile;	//file for writing out image xy coords to binary file
    FILE*                   m_pJointRotFile;    //file for writing out node quaternions

	//Video file writing
	IMFSinkWriter   *m_pSinkWriter;
	DWORD            m_streamIndex;
	LONGLONG         m_rtStart;			//timestamp on video frame

	//Triggering Serial Port
	//static HANDLE	        m_hSerial;

    // Direct2D
	ID2D1Factory*           m_pD2DFactory;
	RGBQUAD*                m_pColorRGBX;
	RGBQUAD*                m_pSwapline;	//for flipping image frame when writing h264 video
	ID2D1Bitmap*            m_pBitmap;

    // Body/hand drawing
    ID2D1HwndRenderTarget*  m_pRenderTarget;
    ID2D1SolidColorBrush*   m_pBrushJointTracked;
    ID2D1SolidColorBrush*   m_pBrushJointInferred;
    ID2D1SolidColorBrush*   m_pBrushBoneTracked;
    ID2D1SolidColorBrush*   m_pBrushBoneInferred;
    ID2D1SolidColorBrush*   m_pBrushHandClosed;
    ID2D1SolidColorBrush*   m_pBrushHandOpen;
    ID2D1SolidColorBrush*   m_pBrushHandLasso;

    /// <summary>
    /// Main processing function
    /// </summary>
    void                    Update();

    /// <summary>
    /// Initializes the default Kinect sensor
    /// </summary>
    /// <returns>S_OK on success, otherwise failure code</returns>
    HRESULT                 InitializeDefaultSensor();
    
    /// <summary>
    /// Handle new body and color data
    /// <param name="nTime">timestamp of body frame</param>
    /// <param name="nBodyCount">body data count</param>
    /// <param name="ppBodies">body data in frame</param>
	/// <param name="nTime">timestamp of color frame</param>
	/// <param name="pBuffer">pointer to frame data</param>
	/// <param name="nWidth">width (in pixels) of input image data</param>
	/// <param name="nHeight">height (in pixels) of input image data</param>
    /// </summary>
	void                    ProcessKinect2Acq(INT64 nTime, int nBodyCount, IBody** ppBodies, RGBQUAD* pBuffer, int nWidth, int nHeight);

    /// <summary>
    /// Set the status bar message
    /// </summary>
    /// <param name="szMessage">message to display</param>
    /// <param name="nShowTimeMsec">time in milliseconds for which to ignore future status messages</param>
    /// <param name="bForce">force status update</param>
    bool                    SetStatusMessage(WCHAR* szMessage, DWORD nShowTimeMsec, bool bForce);

    /// <summary>
    /// Ensure necessary Direct2d resources are created
    /// </summary>
    /// <returns>S_OK if successful, otherwise an error code</returns>
    HRESULT EnsureDirect2DResources();

    /// <summary>
    /// Dispose Direct2d resources 
    /// </summary>
    void DiscardDirect2DResources();

    /// <summary>
    /// Converts a body point to screen space
    /// </summary>
    /// <param name="bodyPoint">body point to tranform</param>
    /// <returns>point in screen-space</returns>
    D2D1_POINT_2F           BodyToScreen(const CameraSpacePoint& bodyPoint);

    /// <summary>
    /// Draws a body 
    /// </summary>
    /// <param name="pJoints">joint data</param>
    /// <param name="pJointPoints">joint positions converted to screen space</param>
    void                    DrawBody(const Joint* pJoints, const D2D1_POINT_2F* pJointPoints);

    /// <summary>
    /// Draws a hand symbol if the hand is tracked: red circle = closed, green circle = opened; blue circle = lasso
    /// </summary>
    /// <param name="handState">state of the hand</param>
    /// <param name="handPosition">position of the hand</param>
    void                    DrawHand(HandState handState, const D2D1_POINT_2F& handPosition);

    /// <summary>
    /// Draws one bone of a body (joint to joint)
    /// </summary>
    /// <param name="pJoints">joint data</param>
    /// <param name="pJointPoints">joint positions converted to screen space</param>
    /// <param name="pJointPoints">joint positions converted to screen space</param>
    /// <param name="joint0">one joint of the bone to draw</param>
    /// <param name="joint1">other joint of the bone to draw</param>
    void                    DrawBone(const Joint* pJoints, const D2D1_POINT_2F* pJointPoints, JointType joint0, JointType joint1);

	/// <summary>
	/// Base name of the data files to be stored upon pressing "Record".
	/// </summary>
	/// <param name="lpszFilePath">string buffer that will receive File name.</param>
	/// <param name="nFilePathSize">number of characters in lpszFilePath string buffer.</param>
	/// <returns>
	/// S_OK on success, otherwise failure code.
	/// </returns>
	HRESULT                 GetTrackingFileName(LPWSTR lpszFilePath, UINT nFilePathSize);

	/// <summary>
	/// Open Files for writing video, skeletal tracking (3d/2d), and joint rotations to file
	/// </summary>
	/// <param name="TrackingDataPath">Base name of files</param>
	/// <returns>indicates success or failure</returns>
	HRESULT                 OpenDataFiles(WCHAR *TrackingDataPath);

	/// <summary>
	/// Write video, skeletal tracking (3d/2d), and joint rotations to file
	/// </summary>
	/// <param name="jointsXYZ">buffer containing joint coords, xyz</param>
	/// <param name="count0">integer count of total joint coords, world xyz coords</param>
	/// <param name="jointsXYZ">buffer containing joint coords, image xy coords</param>
	/// <param name="count1">integer count of total joint coords, xy</param>
    /// <param name="jointRots">buffer containing joint quarternions, xyzw</param>
    /// <param name="count2">integer count of total joint rotations</param>
	/// <returns>indicates success or failure</returns>
	HRESULT                WriteDataFiles(float *jointsXYZ, const int count0, 
                                          float *jointsXY,  const int count1,
                                          float *jointRots, const int count2);

	/// <summary>
	/// Close Data files
	/// </summary>
	/// <returns>indicates success or failure</returns>
	HRESULT                CloseDataFiles();

	/// <summary>
	/// Initialize Sink Writer video writing object
	/// </summary>
	/// <param name="IMFSinkWriter">A pointer to an uninitialize Sink Writer object</param>
	/// <param name="DWORD">A stream index for the video writer</param>
	/// <param name="LPCWSTR">video fiel name to be recorded</param>
	/// <returns>indicates success or failure</returns>
	HRESULT InitializeSinkWriter(IMFSinkWriter **ppWriter, DWORD *pStreamIndex, LPCWSTR videoFilename);

	/// <summary>
	/// Write a video frame using a Sink Writer
	/// </summary>
	/// <param name="IMFSinkWriter">A Sink Writer object</param>
	/// <param name="DWORD">A stream index for the video writer</param>
	/// <param name="DWORD">Time stamp for the frame</param>
	/// <returns>indicates success or failure</returns>
	HRESULT WriteFrame(IMFSinkWriter *pWriter, DWORD streamIndex, const LONGLONG& rtStart);

	/// <summary>
	/// Flip RGBQUAD image frame from bottom-up to top-down
	/// </summary>
	/// <param name="RGBQUAD">An image of RGBQUAD pixels</param>
	/// <returns>flipped image</returns>
	void FlipFrame(RGBQUAD *Frame);

	/// <summary>
	/// Send trigger Pulse via Serial Port DTR pin on separate thread
	/// </summary>
	/// <returns>indicates success or failure</returns>
	//static void                TriggerSerialPort( void * );

	/// <summary>
	/// Get the name of the file where screenshot will be stored.
	/// </summary>
	/// <param name="lpszFilePath">string buffer that will receive screenshot file name.</param>
	/// <param name="nFilePathSize">number of characters in lpszFilePath string buffer.</param>
	/// <returns>
	/// S_OK on success, otherwise failure code.
	/// </returns>
	HRESULT                 GetScreenshotFileName(LPWSTR lpszFilePath, UINT nFilePathSize);

	/// <summary>
	/// Save passed in image data to disk as a bitmap
	/// </summary>
	/// <param name="pBitmapBits">image data to save</param>
	/// <param name="lWidth">width (in pixels) of input image data</param>
	/// <param name="lHeight">height (in pixels) of input image data</param>
	/// <param name="wBitsPerPixel">bits per pixel of image data</param>
	/// <param name="lpszFilePath">full file path to output bitmap to</param>
	/// <returns>indicates success or failure</returns>
	HRESULT                 SaveBitmapToFile(BYTE* pBitmapBits, LONG lWidth, LONG lHeight, WORD wBitsPerPixel, LPCWSTR lpszFilePath);
};

