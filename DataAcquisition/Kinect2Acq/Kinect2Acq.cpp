//------------------------------------------------------------------------------
// <copyright file="Kinect2Acq.cpp" company="Microsoft">
//     The Kinect for Windows APIs used here are preliminary and subject to change
//     Copyright (c) Microsoft Corporation.  All rights reserved.
// </copyright>
//------------------------------------------------------------------------------

#include "stdafx.h"
#include <strsafe.h>
#include "resource.h"
#include "Kinect2Acq.h"
#include <process.h>

static const float c_JointThickness = 3.0f;
static const float c_TrackedBoneThickness = 6.0f;
static const float c_InferredBoneThickness = 1.0f;
static const float c_HandSize = 30.0f;

double CKinect2Acq::m_fFreq = 0.0;

//Global for storing the user input file name
wchar_t* FileBaseName;

/*
//Define static members
HANDLE CKinect2Acq::m_hSerial = CreateFile(L"COM1",
	GENERIC_READ | GENERIC_WRITE,
	0,
	0,
	OPEN_EXISTING,
	FILE_ATTRIBUTE_NORMAL,
	0);

float CKinect2Acq::m_PulseLength = 1.0 / 60.0f;
*/

/// <summary>
/// Entry point for the application
/// </summary>
/// <param name="hInstance">handle to the application instance</param>
/// <param name="hPrevInstance">always 0</param>
/// <param name="lpCmdLine">command line arguments</param>
/// <param name="nCmdShow">whether to display minimized, maximized, or normally</param>
/// <returns>status</returns>
int APIENTRY wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPWSTR lpCmdLine, int nCmdShow)
{	

	UNREFERENCED_PARAMETER(hPrevInstance);
    
	LPWSTR *szArglist;
	int nArgs;

	szArglist = CommandLineToArgvW(GetCommandLineW(), &nArgs);
	if (NULL == szArglist)
	{
		MessageBox(NULL, L"CommandLineToArgvW failed", L"Failed", MB_OK);
		return 0;
	}

	if (nArgs < 2)
	{
		//the user did not enter a base name for the file
		FileBaseName = L"KinematicTracking";
	}
	else if (nArgs > 2)
	{
		//The user put in too many arguments
		MessageBox(NULL, L"Only 1 argument is accepted: File base name", L"Too Many Arguments", MB_OK);
		return 0;
	}
	else
	{
		FileBaseName = szArglist[1];
	}

    CKinect2Acq application;
    application.Run(hInstance, nCmdShow);
}

/// <summary>
/// Constructor
/// </summary>
CKinect2Acq::CKinect2Acq() :
    m_hWnd(NULL),
    m_nStartTime(0),
    m_nLastCounter(0),
    m_nFramesSinceUpdate(0),
    m_nNextStatusTime(0),
    m_pKinectSensor(NULL),
    m_pCoordinateMapper(NULL),
	m_pMultiSourceFrameReader(NULL),
	m_bRecordData(false),
	m_pDataFile(NULL),
	m_pImageCoordsFile(NULL),
	m_pJointRotFile(NULL),
	m_pSinkWriter(NULL),
	m_streamIndex(0),
	m_rtStart(0),
	m_bSaveScreenshot(false),
	m_pD2DFactory(NULL),
	m_pColorRGBX(NULL),
	m_pSwapline(NULL),
	m_pBitmap(0),
    m_pRenderTarget(NULL),
    m_pBrushJointTracked(NULL),
    m_pBrushJointInferred(NULL),
    m_pBrushBoneTracked(NULL),
    m_pBrushBoneInferred(NULL),
    m_pBrushHandClosed(NULL),
    m_pBrushHandOpen(NULL),
    m_pBrushHandLasso(NULL)
{

	LARGE_INTEGER qpf = {0};
    if (QueryPerformanceFrequency(&qpf))
    {
        m_fFreq = double(qpf.QuadPart);
    }

	// create heap storage for color pixel data in RGBX format
	m_pColorRGBX = new RGBQUAD[cColorWidth * cColorHeight];
	//create heap storage for swap line used in flipping image from for video writing
	m_pSwapline = new RGBQUAD[cColorWidth*sizeof(RGBQUAD)];

}
  

/// <summary>
/// Destructor
/// </summary>
CKinect2Acq::~CKinect2Acq()
{
    DiscardDirect2DResources();

	// clean up Direct2D renderer
	if (m_pColorRGBX)
	{
		delete[] m_pColorRGBX;
		m_pColorRGBX = NULL;
	}

	// clean up Direct2D renderer
	if (m_pSwapline)
	{
		delete[] m_pSwapline;
		m_pSwapline = NULL;
	}

	//clean up binary files
	if (m_pDataFile)
	{
		CloseDataFiles();
	}

    // clean up Direct2D
    SafeRelease(m_pD2DFactory);

	// done with frame reader
	SafeRelease(m_pMultiSourceFrameReader);

    // done with coordinate mapper
    SafeRelease(m_pCoordinateMapper);

    // close the Kinect Sensor
    if (m_pKinectSensor)
    {
        m_pKinectSensor->Close();
    }

    SafeRelease(m_pKinectSensor);
}

/// <summary>
/// Creates the main window and begins processing
/// </summary>
/// <param name="hInstance">handle to the application instance</param>
/// <param name="nCmdShow">whether to display minimized, maximized, or normally</param>
int CKinect2Acq::Run(HINSTANCE hInstance, int nCmdShow)
{
		
    MSG       msg = {0};
    WNDCLASS  wc;

    // Dialog custom window class
    ZeroMemory(&wc, sizeof(wc));
    wc.style         = CS_HREDRAW | CS_VREDRAW;
    wc.cbWndExtra    = DLGWINDOWEXTRA;
    wc.hCursor       = LoadCursorW(NULL, IDC_ARROW);
    wc.hIcon         = LoadIconW(hInstance, MAKEINTRESOURCE(IDI_APP));
    wc.lpfnWndProc   = DefDlgProcW;
    wc.lpszClassName = L"Kinect2AcqAppDlgWndClass";

    if (!RegisterClassW(&wc))
    {
        return 0;
    }

    // Create main application window
    HWND hWndApp = CreateDialogParamW(
        NULL,
        MAKEINTRESOURCE(IDD_APP),
        NULL,
        (DLGPROC)CKinect2Acq::MessageRouter, 
        reinterpret_cast<LPARAM>(this));

	// Show window
    ShowWindow(hWndApp, nCmdShow);

	/*
	if (m_hSerial == INVALID_HANDLE_VALUE){
		if (GetLastError() == ERROR_FILE_NOT_FOUND){
			return -1;
		}
		return -2;
	}

	DCB dbcSerialParams = { 0 };

	dbcSerialParams.DCBlength = sizeof(dbcSerialParams);

	if (!GetCommState(m_hSerial, &dbcSerialParams)){
		return -3;
	}

	dbcSerialParams.BaudRate = 9600;
	dbcSerialParams.ByteSize = 8;
	dbcSerialParams.StopBits = ONESTOPBIT;
	dbcSerialParams.Parity = NOPARITY;
	dbcSerialParams.fDtrControl = DTR_CONTROL_DISABLE;
	dbcSerialParams.fRtsControl = RTS_CONTROL_DISABLE;

	if (!SetCommState(m_hSerial, &dbcSerialParams)){
		return -4;
	}
	*/

    // Main message loop
    while (WM_QUIT != msg.message)
    {
        Update();

        while (PeekMessageW(&msg, NULL, 0, 0, PM_REMOVE))
        {
            // If a dialog message will be taken care of by the dialog proc
            if (hWndApp && IsDialogMessageW(hWndApp, &msg))
            {
                continue;
            }

            TranslateMessage(&msg);
            DispatchMessageW(&msg);
        }
    }

    return static_cast<int>(msg.wParam);
}

/// <summary>
/// Main processing function
/// </summary>
void CKinect2Acq::Update()
{
	
	if (!m_pMultiSourceFrameReader)
	{
		SetStatusMessage(L"Failed to initialize MultiSourceFrameReader.", 10000, true);
		return;
	}

	//Query Frames for both skeletal tracking and color image
	IMultiSourceFrame* pMultiSourceFrame = NULL;
    IBodyFrame* pBodyFrame = NULL;
	IColorFrame* pColorFrame = NULL;

	HRESULT hr = m_pMultiSourceFrameReader->AcquireLatestFrame(&pMultiSourceFrame);

	//Get references to individual sensor streams
	if (SUCCEEDED(hr))
    {
		IBodyFrameReference* pBodyFrameReference = NULL;

		hr = pMultiSourceFrame->get_BodyFrameReference(&pBodyFrameReference);
		if (SUCCEEDED(hr))
		{
			hr = pBodyFrameReference->AcquireFrame(&pBodyFrame);
		}

		SafeRelease(pBodyFrameReference);
	}

	if (SUCCEEDED(hr))
	{
		IColorFrameReference* pColorFrameReference = NULL;

		hr = pMultiSourceFrame->get_ColorFrameReference(&pColorFrameReference);
		if (SUCCEEDED(hr))
		{
			hr = pColorFrameReference->AcquireFrame(&pColorFrame);
		}

		SafeRelease(pColorFrameReference);
	}

	if (SUCCEEDED(hr))
	{
        INT64 nTime = 0;
		//Initialize Image frame description
		IFrameDescription* pFrameDescription = NULL;
		int nWidth = 0;
		int nHeight = 0;
		ColorImageFormat imageFormat = ColorImageFormat_None;
		UINT nBufferSize = 0;
		RGBQUAD *pBuffer = NULL;

		// get body frame data
        hr = pBodyFrame->get_RelativeTime(&nTime);
		
        IBody* ppBodies[BODY_COUNT] = {0};

		//Body Specific Checks
        if (SUCCEEDED(hr))
        {
            hr = pBodyFrame->GetAndRefreshBodyData(_countof(ppBodies), ppBodies);
        }

		//Color Image Specific Checks
		if (SUCCEEDED(hr))
		{
			hr = pColorFrame->get_FrameDescription(&pFrameDescription);
		}

		if (SUCCEEDED(hr))
		{
			hr = pFrameDescription->get_Width(&nWidth);
		}

		if (SUCCEEDED(hr))
		{
			hr = pFrameDescription->get_Height(&nHeight);
		}

		if (SUCCEEDED(hr))
		{
			hr = pColorFrame->get_RawColorImageFormat(&imageFormat);
		}

		if (SUCCEEDED(hr))
		{
			if (imageFormat == ColorImageFormat_Bgra)
			{
				hr = pColorFrame->AccessRawUnderlyingBuffer(&nBufferSize, reinterpret_cast<BYTE**>(&pBuffer));
			}
			else if (m_pColorRGBX)
			{
				pBuffer = m_pColorRGBX;
				nBufferSize = cColorWidth * cColorHeight * sizeof(RGBQUAD);
				hr = pColorFrame->CopyConvertedFrameDataToArray(nBufferSize, reinterpret_cast<BYTE*>(pBuffer), ColorImageFormat_Bgra);
			}
			else
			{
				hr = E_FAIL;
			}
		}

		if (SUCCEEDED(hr))
        {
			ProcessKinect2Acq(nTime, BODY_COUNT, ppBodies, pBuffer, nWidth, nHeight);
        }

        for (int i = 0; i < _countof(ppBodies); ++i)
        {
            SafeRelease(ppBodies[i]);
        }

		SafeRelease(pFrameDescription);
    }

    SafeRelease(pBodyFrame);
	SafeRelease(pColorFrame);
	SafeRelease(pMultiSourceFrame);
}

/// <summary>
/// Handles window messages, passes most to the class instance to handle
/// </summary>
/// <param name="hWnd">window message is for</param>
/// <param name="uMsg">message</param>
/// <param name="wParam">message data</param>
/// <param name="lParam">additional message data</param>
/// <returns>result of message processing</returns>
LRESULT CALLBACK CKinect2Acq::MessageRouter(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    CKinect2Acq* pThis = NULL;
    
    if (WM_INITDIALOG == uMsg)
    {
        pThis = reinterpret_cast<CKinect2Acq*>(lParam);
        SetWindowLongPtr(hWnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(pThis));
    }
    else
    {
        pThis = reinterpret_cast<CKinect2Acq*>(::GetWindowLongPtr(hWnd, GWLP_USERDATA));
    }

    if (pThis)
    {
        return pThis->DlgProc(hWnd, uMsg, wParam, lParam);
    }

    return 0;
}

/// <summary>
/// Handle windows messages for the class instance
/// </summary>
/// <param name="hWnd">window message is for</param>
/// <param name="uMsg">message</param>
/// <param name="wParam">message data</param>
/// <param name="lParam">additional message data</param>
/// <returns>result of message processing</returns>
LRESULT CALLBACK CKinect2Acq::DlgProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(wParam);
    UNREFERENCED_PARAMETER(lParam);
	
    switch (message)
    {
        case WM_INITDIALOG:
        {
            // Bind application window handle
            m_hWnd = hWnd;

            // Init Direct2D
            D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, &m_pD2DFactory);
			
			m_sourceStride = cColorWidth*sizeof(RGBQUAD);

            // Get and initialize the default Kinect sensor
            InitializeDefaultSensor();
        }
        break;

        // If the titlebar X is clicked, destroy app
        case WM_CLOSE:
            DestroyWindow(hWnd);
            break;

        case WM_DESTROY:
            // Quit the main message pump
            PostQuitMessage(0);
            break;

			// Handle button press
		case WM_COMMAND:
			// If it was for the Record control and a button clicked event, start recording skeletal data next frame 
			if (IDC_BUTTON_RECORD == LOWORD(wParam) && BN_CLICKED == HIWORD(wParam))
			{
				m_bRecordData = true;
			}
			// If it was for the Stop control and a button clicked event, stop recording skeletal data next frame 
			if (IDC_BUTTON_STOP == LOWORD(wParam) && BN_CLICKED == HIWORD(wParam))
			{
				m_bRecordData = false;
			}
			break;
    }

    return FALSE;
}

/// <summary>
/// Initializes the default Kinect sensor
/// </summary>
/// <returns>indicates success or failure</returns>
HRESULT CKinect2Acq::InitializeDefaultSensor()
{
    HRESULT hr;

    hr = GetDefaultKinectSensor(&m_pKinectSensor);
    if (FAILED(hr))
    {
        return hr;
    }

	if (m_pKinectSensor)
	{
		// Initialize the Kinect and get coordinate mapper and the frame reader

		if (SUCCEEDED(hr))
		{
			hr = m_pKinectSensor->get_CoordinateMapper(&m_pCoordinateMapper);
		}

		hr = m_pKinectSensor->Open();

		if (SUCCEEDED(hr))
		{
			hr = m_pKinectSensor->OpenMultiSourceFrameReader(
				FrameSourceTypes::FrameSourceTypes_Depth | FrameSourceTypes::FrameSourceTypes_Color | FrameSourceTypes::FrameSourceTypes_Body,
				&m_pMultiSourceFrameReader);
		}
	}
    if (!m_pKinectSensor || FAILED(hr))
    {
        SetStatusMessage(L"No ready Kinect found!", 10000, true);
        return E_FAIL;
    }

    return hr;
}

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
void CKinect2Acq::ProcessKinect2Acq(INT64 nTime, int nBodyCount, IBody** ppBodies, RGBQUAD* pBuffer, int nWidth, int nHeight)
{
    if (m_hWnd)
    {

		HRESULT hr = EnsureDirect2DResources();

		//DRAW RGB IMAGE FIRST
		// Make sure we've received valid data
		if (SUCCEEDED(hr) && m_pRenderTarget && m_pCoordinateMapper &&
			pBuffer && (nWidth == cColorWidth) && (nHeight == cColorHeight))
		{
			
			// Draw the data with Direct2D
			// Copy the image that was passed in into the direct2d bitmap
			hr = m_pBitmap->CopyFromMemory(NULL, reinterpret_cast<BYTE*>(pBuffer), m_sourceStride);

			if (FAILED(hr))
			{
				SetStatusMessage(L"Failed to Copy Image Buffer to Bitmap.", 10000, true);
				return;
			}

			//Initialize Drawing
			m_pRenderTarget->BeginDraw();
			// Draw the bitmap stretched to the size of the window
			m_pRenderTarget->DrawBitmap(m_pBitmap);
			
			//DRAW SKELETON
			for (int i = 0; i < nBodyCount; ++i)
			{
				IBody* pBody = ppBodies[i];
				if (pBody)
				{

					BOOLEAN bTracked = false;
					hr = pBody->get_IsTracked(&bTracked);

					if (SUCCEEDED(hr) && bTracked)
					{
						/*
						//Trigger Serial Port DTR pin
						if (m_pDataFile)
						{
								hr = _beginthread(TriggerSerialPort, 0, (void *)m_hSerial);

								if (FAILED(hr))
								{
									SetStatusMessage(L"Failed to Trigger Serial Port.", 10000, true);
									return;
								}
						}
						*/

						//Access Joint Data
						Joint joints[JointType_Count];
						D2D1_POINT_2F jointPoints[JointType_Count];
						JointOrientation jointOs[JointType_Count];

						/*Buffers for writing out:
						1) xyz world coords of each node
						2) xy image plane coords of each node
						3) 4 component quaternion detailing the rotation angle of each node
						*/
						float jointBuffer[3*JointType_Count];
						float jointXY[2*JointType_Count];
						float jointRots[4*JointType_Count];


						HandState leftHandState = HandState_Unknown;
						HandState rightHandState = HandState_Unknown;

						pBody->get_HandLeftState(&leftHandState);
						pBody->get_HandRightState(&rightHandState);
						
						hr = pBody->GetJoints(_countof(joints), joints);
						hr = pBody->GetJointOrientations(_countof(jointOs), jointOs);
						
						if (SUCCEEDED(hr))
						{
							for (int j = 0; j < _countof(joints); ++j)
							{
								//Transform this joint into Screen space
								jointPoints[j] = BodyToScreen(joints[j].Position);
								//Write this nodes data to it's appropriate array locations
								//world frame xyz
								jointBuffer[3*joints[j].JointType]   = joints[j].Position.X;
								jointBuffer[3*joints[j].JointType+1] = joints[j].Position.Y;
								jointBuffer[3*joints[j].JointType+2] = joints[j].Position.Z;
								//image plane xy
								jointXY[2*joints[j].JointType]       = jointPoints[j].x;
								jointXY[2*joints[j].JointType+1]     = jointPoints[j].y;
								//node rotation quaternion
								jointRots[4*jointOs[j].JointType]	 = jointOs[j].Orientation.x;
								jointRots[4*jointOs[j].JointType+1]	 = jointOs[j].Orientation.y;
								jointRots[4*jointOs[j].JointType+2]	 = jointOs[j].Orientation.z;
								jointRots[4*jointOs[j].JointType+3]	 = jointOs[j].Orientation.w;
							}

							DrawBody(joints, jointPoints);

							DrawHand(leftHandState, jointPoints[JointType_HandLeft]);
							DrawHand(rightHandState, jointPoints[JointType_HandRight]);

							//If files exist, write data to files
							if (m_pDataFile){
								hr = WriteDataFiles(jointBuffer, 3*JointType_Count, 
													jointXY,     2*JointType_Count,
													jointRots,   4*JointType_Count);
							}

							if (FAILED(hr))
							{
								SetStatusMessage(L"Failed to Write to binary file.", 10000, true);
								return;
							}
						}
					}
				}
			}
			
			hr = m_pRenderTarget->EndDraw();

			// Device lost, need to recreate the render target
			// We'll dispose it now and retry drawing
			if (D2DERR_RECREATE_TARGET == hr)
			{
				hr = S_OK;
				DiscardDirect2DResources();
			}
		}

        if (!m_nStartTime)
        {
            m_nStartTime = nTime;
        }

        double fps = 0.0;

        LARGE_INTEGER qpcNow = {0};
        if (m_fFreq)
        {
            if (QueryPerformanceCounter(&qpcNow))
            {
                if (m_nLastCounter)
                {
                    m_nFramesSinceUpdate++;
                    fps = m_fFreq * m_nFramesSinceUpdate / double(qpcNow.QuadPart - m_nLastCounter);
                }
            }
        }

        WCHAR szStatusMessage[64];
        StringCchPrintf(szStatusMessage, _countof(szStatusMessage), L" FPS = %0.2f    Time = %I64d", fps, (nTime - m_nStartTime));

        if (SetStatusMessage(szStatusMessage, 1000, false))
        {
            m_nLastCounter = qpcNow.QuadPart;
            m_nFramesSinceUpdate = 0;
        }

		//Start Recording Skeletal Data
		if (m_bRecordData && !m_pDataFile)
		{
			WCHAR szTrackingDataPath[MAX_PATH];

			// Retrieve the path to My Videos
			GetTrackingFileName(szTrackingDataPath, _countof(szTrackingDataPath));

			// Open binary file for writing skeletal tracking data
			HRESULT hr = OpenDataFiles(szTrackingDataPath);
			
			WCHAR szStatusMessage[64 + MAX_PATH];
			if (SUCCEEDED(hr))
			{
				// Set the status bar to show where the screenshot was saved
				StringCchPrintf(szStatusMessage, _countof(szStatusMessage), L"Recording to %s", szTrackingDataPath);
			}
			else
			{
				StringCchPrintf(szStatusMessage, _countof(szStatusMessage), L"Failed to open all data files for %s", szTrackingDataPath);
			}

			SetStatusMessage(szStatusMessage, 3000, true);

		}

		//Stop Recording Skeletal Data
		if (!m_bRecordData && m_pDataFile)
		{
			
			// Open binary file for writing skeletal tracking data
			HRESULT hr = CloseDataFiles();

			WCHAR szStatusMessage[64 + MAX_PATH];
			if (SUCCEEDED(hr))
			{
				// Set the status bar to show where the screenshot was saved
				StringCchPrintf(szStatusMessage, _countof(szStatusMessage), L"Stopped Recording");
			}
			else
			{
				StringCchPrintf(szStatusMessage, _countof(szStatusMessage), L"Failed to close data files!");
			}

			SetStatusMessage(szStatusMessage, 3000, true);

			//toggle binary file state
			m_pDataFile = NULL;
		}

		//Save Screenshot routines
		if (m_bSaveScreenshot)
		{
			WCHAR szScreenshotPath[MAX_PATH];

			// Retrieve the path to My Photos
			GetScreenshotFileName(szScreenshotPath, _countof(szScreenshotPath));

			// Write out the bitmap to disk
			HRESULT hr = SaveBitmapToFile(reinterpret_cast<BYTE*>(pBuffer), nWidth, nHeight, sizeof(RGBQUAD) * 8, szScreenshotPath);

			WCHAR szStatusMessage[64 + MAX_PATH];
			if (SUCCEEDED(hr))
			{
				// Set the status bar to show where the screenshot was saved
				StringCchPrintf(szStatusMessage, _countof(szStatusMessage), L"Screenshot saved to %s", szScreenshotPath);
			}
			else
			{
				StringCchPrintf(szStatusMessage, _countof(szStatusMessage), L"Failed to write screenshot to %s", szScreenshotPath);
			}

			SetStatusMessage(szStatusMessage, 5000, true);

			// toggle off so we don't save a screenshot again next frame
			m_bSaveScreenshot = false;
		}
    }
}

/// <summary>
/// Set the status bar message
/// </summary>
/// <param name="szMessage">message to display</param>
/// <param name="showTimeMsec">time in milliseconds to ignore future status messages</param>
/// <param name="bForce">force status update</param>
bool CKinect2Acq::SetStatusMessage(WCHAR* szMessage, DWORD nShowTimeMsec, bool bForce)
{
    DWORD now = GetTickCount();

    if (m_hWnd && (bForce || (m_nNextStatusTime <= now)))
    {
        SetDlgItemText(m_hWnd, IDC_STATUS, szMessage);
        m_nNextStatusTime = now + nShowTimeMsec;

        return true;
    }

    return false;
}

/// <summary>
/// Ensure necessary Direct2d resources are created
/// </summary>
/// <returns>S_OK if successful, otherwise an error code</returns>
HRESULT CKinect2Acq::EnsureDirect2DResources()
{
    HRESULT hr = S_OK;

    if (m_pD2DFactory && !m_pRenderTarget)
    {
		D2D1_SIZE_U size = D2D1::SizeU(cColorWidth, cColorHeight);

		D2D1_RENDER_TARGET_PROPERTIES rtProps = D2D1::RenderTargetProperties();
		rtProps.pixelFormat = D2D1::PixelFormat(DXGI_FORMAT_B8G8R8A8_UNORM, D2D1_ALPHA_MODE_IGNORE);
		rtProps.usage = D2D1_RENDER_TARGET_USAGE_GDI_COMPATIBLE;

		// Create a hWnd render target, in order to render to the window set in initialize
		hr = m_pD2DFactory->CreateHwndRenderTarget(
			rtProps,
			D2D1::HwndRenderTargetProperties(GetDlgItem(m_hWnd, IDC_VIDEOVIEW), size),
			&m_pRenderTarget
			);

        if (FAILED(hr))
        {
            SetStatusMessage(L"Couldn't create Direct2D render target!", 10000, true);
            return hr;
        }

		// Create a bitmap that we can copy image data into and then render to the target
		hr = m_pRenderTarget->CreateBitmap(
			size,
			D2D1::BitmapProperties(D2D1::PixelFormat(DXGI_FORMAT_B8G8R8A8_UNORM, D2D1_ALPHA_MODE_IGNORE)),
			&m_pBitmap
			);

        // light green
        m_pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(0.27f, 0.75f, 0.27f), &m_pBrushJointTracked);

        m_pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::Yellow, 1.0f), &m_pBrushJointInferred);
        m_pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::Green, 1.0f), &m_pBrushBoneTracked);
        m_pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::Gray, 1.0f), &m_pBrushBoneInferred);

        m_pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::Red, 0.5f), &m_pBrushHandClosed);
        m_pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::Green, 0.5f), &m_pBrushHandOpen);
        m_pRenderTarget->CreateSolidColorBrush(D2D1::ColorF(D2D1::ColorF::Blue, 0.5f), &m_pBrushHandLasso);
		

		if (FAILED(hr))
		{
			SafeRelease(m_pRenderTarget);
			return hr;
		}
    }

    return hr;
}

/// <summary>
/// Dispose Direct2d resources 
/// </summary>
void CKinect2Acq::DiscardDirect2DResources()
{
    SafeRelease(m_pRenderTarget);

    SafeRelease(m_pBrushJointTracked);
    SafeRelease(m_pBrushJointInferred);
    SafeRelease(m_pBrushBoneTracked);
    SafeRelease(m_pBrushBoneInferred);

    SafeRelease(m_pBrushHandClosed);
    SafeRelease(m_pBrushHandOpen);
    SafeRelease(m_pBrushHandLasso);

	SafeRelease(m_pBitmap);
}

/// <summary>
/// Converts a body point to screen space
/// </summary>
/// <param name="bodyPoint">body point to tranform</param>
/// <returns>point in screen-space</returns>
D2D1_POINT_2F CKinect2Acq::BodyToScreen(const CameraSpacePoint& bodyPoint)
{
    // Calculate the body's position on the screen
    ColorSpacePoint colorPoint = {0};
    m_pCoordinateMapper->MapCameraPointToColorSpace(bodyPoint, &colorPoint);

	float screenPointX = static_cast<float>(colorPoint.X);
    float screenPointY = static_cast<float>(colorPoint.Y);

    return D2D1::Point2F(screenPointX, screenPointY);
}

/// <summary>
/// Draws a body 
/// </summary>
/// <param name="pJoints">joint data</param>
/// <param name="pJointPoints">joint positions converted to screen space</param>
void CKinect2Acq::DrawBody(const Joint* pJoints, const D2D1_POINT_2F* pJointPoints)
{
    // Draw the bones

    // Torso
    DrawBone(pJoints, pJointPoints, JointType_Head, JointType_Neck);
    DrawBone(pJoints, pJointPoints, JointType_Neck, JointType_SpineShoulder);
    DrawBone(pJoints, pJointPoints, JointType_SpineShoulder, JointType_SpineMid);
    DrawBone(pJoints, pJointPoints, JointType_SpineMid, JointType_SpineBase);
    DrawBone(pJoints, pJointPoints, JointType_SpineShoulder, JointType_ShoulderRight);
    DrawBone(pJoints, pJointPoints, JointType_SpineShoulder, JointType_ShoulderLeft);
    DrawBone(pJoints, pJointPoints, JointType_SpineBase, JointType_HipRight);
    DrawBone(pJoints, pJointPoints, JointType_SpineBase, JointType_HipLeft);
    
    // Right Arm    
    DrawBone(pJoints, pJointPoints, JointType_ShoulderRight, JointType_ElbowRight);
    DrawBone(pJoints, pJointPoints, JointType_ElbowRight, JointType_WristRight);
    DrawBone(pJoints, pJointPoints, JointType_WristRight, JointType_HandRight);
    DrawBone(pJoints, pJointPoints, JointType_HandRight, JointType_HandTipRight);
    DrawBone(pJoints, pJointPoints, JointType_WristRight, JointType_ThumbRight);

    // Left Arm
    DrawBone(pJoints, pJointPoints, JointType_ShoulderLeft, JointType_ElbowLeft);
    DrawBone(pJoints, pJointPoints, JointType_ElbowLeft, JointType_WristLeft);
    DrawBone(pJoints, pJointPoints, JointType_WristLeft, JointType_HandLeft);
    DrawBone(pJoints, pJointPoints, JointType_HandLeft, JointType_HandTipLeft);
    DrawBone(pJoints, pJointPoints, JointType_WristLeft, JointType_ThumbLeft);

    // Right Leg
    DrawBone(pJoints, pJointPoints, JointType_HipRight, JointType_KneeRight);
    DrawBone(pJoints, pJointPoints, JointType_KneeRight, JointType_AnkleRight);
    DrawBone(pJoints, pJointPoints, JointType_AnkleRight, JointType_FootRight);

    // Left Leg
    DrawBone(pJoints, pJointPoints, JointType_HipLeft, JointType_KneeLeft);
    DrawBone(pJoints, pJointPoints, JointType_KneeLeft, JointType_AnkleLeft);
    DrawBone(pJoints, pJointPoints, JointType_AnkleLeft, JointType_FootLeft);

    // Draw the joints
    for (int i = 0; i < JointType_Count; ++i)
    {
        D2D1_ELLIPSE ellipse = D2D1::Ellipse(pJointPoints[i], c_JointThickness, c_JointThickness);

        if (pJoints[i].TrackingState == TrackingState_Inferred)
        {
            m_pRenderTarget->FillEllipse(ellipse, m_pBrushJointInferred);
        }
        else if (pJoints[i].TrackingState == TrackingState_Tracked)
        {
            m_pRenderTarget->FillEllipse(ellipse, m_pBrushJointTracked);
        }
    }
}

/// <summary>
/// Draws one bone of a body (joint to joint)
/// </summary>
/// <param name="pJoints">joint data</param>
/// <param name="pJointPoints">joint positions converted to screen space</param>
/// <param name="pJointPoints">joint positions converted to screen space</param>
/// <param name="joint0">one joint of the bone to draw</param>
/// <param name="joint1">other joint of the bone to draw</param>
void CKinect2Acq::DrawBone(const Joint* pJoints, const D2D1_POINT_2F* pJointPoints, JointType joint0, JointType joint1)
{
    TrackingState joint0State = pJoints[joint0].TrackingState;
    TrackingState joint1State = pJoints[joint1].TrackingState;

    // If we can't find either of these joints, exit
    if ((joint0State == TrackingState_NotTracked) || (joint1State == TrackingState_NotTracked))
    {
        return;
    }

    // Don't draw if both points are inferred
    if ((joint0State == TrackingState_Inferred) && (joint1State == TrackingState_Inferred))
    {
        return;
    }

    // We assume all drawn bones are inferred unless BOTH joints are tracked
    if ((joint0State == TrackingState_Tracked) && (joint1State == TrackingState_Tracked))
    {
        m_pRenderTarget->DrawLine(pJointPoints[joint0], pJointPoints[joint1], m_pBrushBoneTracked, c_TrackedBoneThickness);
    }
    else
    {
        m_pRenderTarget->DrawLine(pJointPoints[joint0], pJointPoints[joint1], m_pBrushBoneInferred, c_InferredBoneThickness);
    }
}

/// <summary>
/// Draws a hand symbol if the hand is tracked: red circle = closed, green circle = opened; blue circle = lasso
/// </summary>
/// <param name="handState">state of the hand</param>
/// <param name="handPosition">position of the hand</param>
void CKinect2Acq::DrawHand(HandState handState, const D2D1_POINT_2F& handPosition)
{
    D2D1_ELLIPSE ellipse = D2D1::Ellipse(handPosition, c_HandSize, c_HandSize);

    switch (handState)
    {
        case HandState_Closed:
            m_pRenderTarget->FillEllipse(ellipse, m_pBrushHandClosed);
            break;

        case HandState_Open:
            m_pRenderTarget->FillEllipse(ellipse, m_pBrushHandOpen);
            break;

        case HandState_Lasso:
            m_pRenderTarget->FillEllipse(ellipse, m_pBrushHandLasso);
            break;
    }
}

/// <summary>
/// Base name of the data files to be stored upon pressing "Record".
/// </summary>
/// <param name="lpszFilePath">string buffer that will receive File name.</param>
/// <param name="nFilePathSize">number of characters in lpszFilePath string buffer.</param>
/// <returns>
/// S_OK on success, otherwise failure code.
/// </returns>
HRESULT CKinect2Acq::GetTrackingFileName(LPWSTR lpszFilePath, UINT nFilePathSize)
{
	WCHAR* pszKnownPath = NULL;
	HRESULT hr = SHGetKnownFolderPath(FOLDERID_Videos, 0, NULL, &pszKnownPath);

	if (SUCCEEDED(hr))
	{
		//Date string formatting example
		//2015-02-25_141307

		// Get the Date
		WCHAR szDateString[MAX_PATH];
		GetDateFormatEx(NULL, 0, NULL, L"yyyy'-'MM'-'dd", szDateString, _countof(szDateString),NULL);

		// Get the time
		WCHAR szTimeString[MAX_PATH];
		GetTimeFormatEx(NULL, 0, NULL, L"hhmmss", szTimeString, _countof(szTimeString));

		// File name will be FileBaseName_yyyyMMdd-hhmmssX, X= {_3D.bin _2D.bin _Rot.bin .mp4}
		StringCchPrintfW(lpszFilePath, nFilePathSize, L"%s\\%s_%s-%s", pszKnownPath, FileBaseName, szDateString, szTimeString);
	}

	if (pszKnownPath)
	{
		CoTaskMemFree(pszKnownPath);
	}

	return hr;
}

/// <summary>
/// Open Files for writing video, skeletal tracking (3d/2d), and joint rotations to file
/// </summary>
/// <param name="TrackingDataPath">Base name of files</param>
/// <returns>indicates success or failure</returns>
HRESULT CKinect2Acq::OpenDataFiles(WCHAR *szTrackingDataPath)
{
	//Allocate memory for final name of files
	WCHAR szTracking3DPath[MAX_PATH];
	WCHAR szTracking2DPath[MAX_PATH];
	WCHAR szTrackingRotPath[MAX_PATH];
	WCHAR szTrackingVideoPath[MAX_PATH];

	//Open world frame binary file for writing
	StringCchPrintfW(szTracking3DPath, _countof(szTracking3DPath), L"%s_3D.bin", szTrackingDataPath);
	m_pDataFile = _wfopen(szTracking3DPath, L"wb");

	//Return if error opening file
	if (NULL == m_pDataFile)
	{
		return E_ACCESSDENIED;
	}

	//Open image frame binary file for writing
	StringCchPrintfW(szTracking2DPath, _countof(szTracking2DPath), L"%s_2D.bin", szTrackingDataPath);
	m_pImageCoordsFile = _wfopen(szTracking2DPath, L"wb");

	//Return if error opening file
	if (NULL == m_pImageCoordsFile)
	{
		return E_ACCESSDENIED;
	}

	//Open image frame binary file for writing
	StringCchPrintfW(szTrackingRotPath, _countof(szTrackingRotPath), L"%s_Rot.bin", szTrackingDataPath);
	m_pJointRotFile = _wfopen(szTrackingRotPath, L"wb");

	//Return if error opening file
	if (NULL == m_pJointRotFile)
	{
		return E_ACCESSDENIED;
	}

	//Open Sink Writer for writing out video data
	//Initialize Sink Writer
	HRESULT hr = CoInitializeEx(NULL, COINIT_APARTMENTTHREADED);
	if (SUCCEEDED(hr))
	{
		hr = MFStartup(MF_VERSION);
		if (SUCCEEDED(hr))
		{
			StringCchPrintfW(szTrackingVideoPath, _countof(szTrackingVideoPath), L"%s.mp4", szTrackingDataPath);
			hr = InitializeSinkWriter(&m_pSinkWriter, &m_streamIndex, szTrackingVideoPath);
			//reset video frame timestamp
			m_rtStart = 0;
		}
	}

	if (FAILED(hr))
	{
		return E_FAIL;
	}

	return S_OK;
}

HRESULT CKinect2Acq::WriteDataFiles(float *jointsXYZ, const int count0, 
									float *jointsXY,  const int count1,
									float *jointRots, const int count2)
{
	size_t result;

	//write out camera xyz coords for this frame to binary file
	result = fwrite(jointsXYZ, sizeof(float), count0, m_pDataFile);
	//Return if the requested number of elements weren't written
	if (result != count0){
		return E_FAIL;
	}

	//write out image xy coords for this frame to binary file
	result = fwrite(jointsXY, sizeof(float), count1, m_pImageCoordsFile);
	//Return if the requested number of elements weren't written
	if (result != count1){
		return E_FAIL;
	}

	//write out node quarternions for this frame to binary file
	result = fwrite(jointRots, sizeof(float), count2, m_pJointRotFile);
	//Return if the requested number of elements weren't written
	if (result != count2){
		return E_FAIL;
	}

	//write out video frame
	HRESULT hr = WriteFrame(m_pSinkWriter, m_streamIndex, m_rtStart);
	if (FAILED(hr))
	{
		return E_FAIL;
	}
	//increment video frame timestamp
	m_rtStart += VIDEO_FRAME_DURATION;

	return S_OK;
}

HRESULT CKinect2Acq::CloseDataFiles()
{
	int chk;

	//close worlds frame coords binary file
	chk = fclose(m_pDataFile);
	//Return is the file failed to close
	if (chk != 0){
		return E_FAIL;
	}

	//close image frame coords binary file
	chk = fclose(m_pImageCoordsFile);
	//Return is the file failed to close
	if (chk != 0){
		return E_FAIL;
	}

	//close quarternion binary file
	chk = fclose(m_pJointRotFile);
	//Return is the file failed to close
	if (chk != 0){
		return E_FAIL;
	}

	//Finish writing video file and close Sink Writer
	HRESULT hr = m_pSinkWriter->Finalize();

	if (FAILED(hr))
	{
		return E_FAIL;
	}
	SafeRelease(&m_pSinkWriter);
	MFShutdown();

	CoUninitialize();

	return S_OK;
}

HRESULT CKinect2Acq::InitializeSinkWriter(IMFSinkWriter **ppWriter, DWORD *pStreamIndex, LPCWSTR videoFileName)
{
	*ppWriter = NULL;
	*pStreamIndex = NULL;

	IMFSinkWriter   *pSinkWriter = NULL;
	IMFMediaType    *pMediaTypeOut = NULL;
	IMFMediaType    *pMediaTypeIn = NULL;
	DWORD           streamIndex;

	HRESULT hr = MFCreateSinkWriterFromURL(videoFileName, NULL, NULL, &pSinkWriter);

	// Set the output media type.
	if (SUCCEEDED(hr))
	{
		hr = MFCreateMediaType(&pMediaTypeOut);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeOut->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeOut->SetGUID(MF_MT_SUBTYPE, VIDEO_ENCODING_FORMAT);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeOut->SetUINT32(MF_MT_AVG_BITRATE, VIDEO_BIT_RATE);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeOut->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeOut->SetUINT32(MF_MT_MPEG2_PROFILE, eAVEncH264VProfile_Main);
	}
	if (SUCCEEDED(hr))
	{
		hr = MFSetAttributeSize(pMediaTypeOut, MF_MT_FRAME_SIZE, VIDEO_WIDTH, VIDEO_HEIGHT);
	}
	if (SUCCEEDED(hr))
	{
		hr = MFSetAttributeRatio(pMediaTypeOut, MF_MT_FRAME_RATE, VIDEO_FPS, 1);
	}
	if (SUCCEEDED(hr))
	{
		hr = MFSetAttributeRatio(pMediaTypeOut, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);
	}
	if (SUCCEEDED(hr))
	{
		hr = pSinkWriter->AddStream(pMediaTypeOut, &streamIndex);
	}

	// Set the input media type.
	if (SUCCEEDED(hr))
	{
		hr = MFCreateMediaType(&pMediaTypeIn);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeIn->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Video);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeIn->SetGUID(MF_MT_SUBTYPE, VIDEO_INPUT_FORMAT);
	}
	if (SUCCEEDED(hr))
	{
		hr = pMediaTypeIn->SetUINT32(MF_MT_INTERLACE_MODE, MFVideoInterlace_Progressive);
	}
	if (SUCCEEDED(hr))
	{
		hr = MFSetAttributeSize(pMediaTypeIn, MF_MT_FRAME_SIZE, VIDEO_WIDTH, VIDEO_HEIGHT);
	}
	if (SUCCEEDED(hr))
	{
		hr = MFSetAttributeRatio(pMediaTypeIn, MF_MT_FRAME_RATE, VIDEO_FPS, 1);
	}
	if (SUCCEEDED(hr))
	{
		hr = MFSetAttributeRatio(pMediaTypeIn, MF_MT_PIXEL_ASPECT_RATIO, 1, 1);
	}
	if (SUCCEEDED(hr))
	{
		hr = pSinkWriter->SetInputMediaType(streamIndex, pMediaTypeIn, NULL);
	}

	// Tell the sink writer to start accepting data.
	if (SUCCEEDED(hr))
	{
		hr = pSinkWriter->BeginWriting();
	}

	// Return the pointer to the caller.
	if (SUCCEEDED(hr))
	{
		*ppWriter = pSinkWriter;
		(*ppWriter)->AddRef();
		*pStreamIndex = streamIndex;
	}

	SafeRelease(&pSinkWriter);
	SafeRelease(&pMediaTypeOut);
	SafeRelease(&pMediaTypeIn);

	return hr;
}

/// <summary>
/// Write a video frame using a Sink Writer
/// </summary>
/// <param name="IMFSinkWriter">A Sink Writer object</param>
/// <param name="DWORD">A stream index for the video writer</param>
/// <param name="DWORD">Time stamp for the frame</param>
/// <returns>indicates success or failure</returns>
HRESULT CKinect2Acq::WriteFrame(IMFSinkWriter *pWriter, DWORD streamIndex, const LONGLONG& rtStart)
{
	IMFSample *pSample      = NULL;
	IMFMediaBuffer *pBuffer = NULL;
	
	const LONG cbWidth   = m_sourceStride;
	const DWORD cbBuffer = cbWidth * VIDEO_HEIGHT;

	BYTE *pData = NULL;

	// Create a new memory buffer.
	HRESULT hr = MFCreateMemoryBuffer(cbBuffer, &pBuffer);

	// Lock the buffer and copy the video frame to the buffer.
	if (SUCCEEDED(hr))
	{
		hr = pBuffer->Lock(&pData, NULL, NULL);
	}

	//Flip the image frame prior to writing
	FlipFrame(m_pColorRGBX);

	if (SUCCEEDED(hr))
	{
		hr = MFCopyImage(
			pData,									  // Destination buffer.
			cbWidth,								  // Destination stride.
			reinterpret_cast<BYTE*>(m_pColorRGBX),    // First row in source image.
			cbWidth,							      // Source stride.
			cbWidth,                                  // Image width in bytes.
			VIDEO_HEIGHT                              // Image height in pixels.
			);
	}
	if (pBuffer)
	{
		pBuffer->Unlock();
	}

	// Set the data length of the buffer.
	if (SUCCEEDED(hr))
	{
		hr = pBuffer->SetCurrentLength(cbBuffer);
	}

	// Create a media sample and add the buffer to the sample.
	if (SUCCEEDED(hr))
	{
		hr = MFCreateSample(&pSample);
	}
	if (SUCCEEDED(hr))
	{
		hr = pSample->AddBuffer(pBuffer);
	}

	// Set the time stamp and the duration.
	if (SUCCEEDED(hr))
	{
		hr = pSample->SetSampleTime(rtStart);
	}
	if (SUCCEEDED(hr))
	{
		hr = pSample->SetSampleDuration(VIDEO_FRAME_DURATION);
	}

	// Send the sample to the Sink Writer.
	if (SUCCEEDED(hr))
	{
		hr = pWriter->WriteSample(streamIndex, pSample);
	}

	SafeRelease(&pSample);
	SafeRelease(&pBuffer);
	return hr;
}

/// <summary>
/// Flip RGBQUAD image frame from bottom-up to top-down
/// </summary>
/// <param name="RGBQUAD">An image of RGBQUAD pixels</param>
void CKinect2Acq::FlipFrame(RGBQUAD *Frame)
{
	int nRow;
	for (nRow = 0; nRow < cColorHeight / 2; nRow++)
	{
		memcpy(m_pSwapline, (BYTE*)Frame + nRow * m_sourceStride, m_sourceStride);
		memcpy((BYTE*)Frame + nRow*m_sourceStride, (BYTE*)Frame + (cColorHeight - nRow - 1)*m_sourceStride, m_sourceStride);
		memcpy((BYTE*)Frame + (cColorHeight - nRow - 1) * m_sourceStride, m_pSwapline, m_sourceStride);
	}

}

/*
/// <summary>
/// Send trigger Pulse via Serial Port DTR pin on separate thread
/// </summary>
void CKinect2Acq::TriggerSerialPort(void * portHandle)
{
	LARGE_INTEGER start = { 0 };
	LARGE_INTEGER end   = { 0 };
	double dt			= 0.0f;
	
	QueryPerformanceCounter((LARGE_INTEGER*)&start);
	if (!EscapeCommFunction(portHandle, SETDTR)){
		//cout << "Error in setting DTR" << endl;
		//return -6;
	}
	//let it run for PulseLength
	while (dt < m_PulseLength){
		QueryPerformanceCounter(&end);
		dt = double(end.QuadPart - start.QuadPart) / m_fFreq;
	}
	if (!EscapeCommFunction(portHandle, CLRDTR)){
		//cout << "Error clearing DTR" << endl;
		//return -7;
	}
	// _endthread given to terminate
	_endthread();
};
*/

/// <summary>
/// Get the name of the file where screenshot will be stored.
/// </summary>
/// <param name="lpszFilePath">string buffer that will receive screenshot file name.</param>
/// <param name="nFilePathSize">number of characters in lpszFilePath string buffer.</param>
/// <returns>
/// S_OK on success, otherwise failure code.
/// </returns>
HRESULT CKinect2Acq::GetScreenshotFileName(LPWSTR lpszFilePath, UINT nFilePathSize)
{
	WCHAR* pszKnownPath = NULL;
	HRESULT hr = SHGetKnownFolderPath(FOLDERID_Pictures, 0, NULL, &pszKnownPath);

	if (SUCCEEDED(hr))
	{
		// Get the time
		WCHAR szTimeString[MAX_PATH];
		GetTimeFormatEx(NULL, 0, NULL, L"hh'-'mm'-'ss", szTimeString, _countof(szTimeString));

		// File name will be KinectScreenshotColor-HH-MM-SS.bmp
		StringCchPrintfW(lpszFilePath, nFilePathSize, L"%s\\KinectScreenshot-Color-%s.bmp", pszKnownPath, szTimeString);
	}

	if (pszKnownPath)
	{
		CoTaskMemFree(pszKnownPath);
	}

	return hr;
}

/// <summary>
/// Save passed in image data to disk as a bitmap
/// </summary>
/// <param name="pBitmapBits">image data to save</param>
/// <param name="lWidth">width (in pixels) of input image data</param>
/// <param name="lHeight">height (in pixels) of input image data</param>
/// <param name="wBitsPerPixel">bits per pixel of image data</param>
/// <param name="lpszFilePath">full file path to output bitmap to</param>
/// <returns>indicates success or failure</returns>
HRESULT CKinect2Acq::SaveBitmapToFile(BYTE* pBitmapBits, LONG lWidth, LONG lHeight, WORD wBitsPerPixel, LPCWSTR lpszFilePath)
{
	DWORD dwByteCount = lWidth * lHeight * (wBitsPerPixel / 8);

	BITMAPINFOHEADER bmpInfoHeader = { 0 };

	bmpInfoHeader.biSize = sizeof(BITMAPINFOHEADER);  // Size of the header
	bmpInfoHeader.biBitCount = wBitsPerPixel;             // Bit count
	bmpInfoHeader.biCompression = BI_RGB;                    // Standard RGB, no compression
	bmpInfoHeader.biWidth = lWidth;                    // Width in pixels
	bmpInfoHeader.biHeight = -lHeight;                  // Height in pixels, negative indicates it's stored right-side-up
	bmpInfoHeader.biPlanes = 1;                         // Default
	bmpInfoHeader.biSizeImage = dwByteCount;               // Image size in bytes

	BITMAPFILEHEADER bfh = { 0 };

	bfh.bfType = 0x4D42;                                           // 'M''B', indicates bitmap
	bfh.bfOffBits = bmpInfoHeader.biSize + sizeof(BITMAPFILEHEADER);  // Offset to the start of pixel data
	bfh.bfSize = bfh.bfOffBits + bmpInfoHeader.biSizeImage;        // Size of image + headers

	// Create the file on disk to write to
	HANDLE hFile = CreateFileW(lpszFilePath, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);

	// Return if error opening file
	if (NULL == hFile)
	{
		return E_ACCESSDENIED;
	}

	DWORD dwBytesWritten = 0;

	// Write the bitmap file header
	if (!WriteFile(hFile, &bfh, sizeof(bfh), &dwBytesWritten, NULL))
	{
		CloseHandle(hFile);
		return E_FAIL;
	}

	// Write the bitmap info header
	if (!WriteFile(hFile, &bmpInfoHeader, sizeof(bmpInfoHeader), &dwBytesWritten, NULL))
	{
		CloseHandle(hFile);
		return E_FAIL;
	}

	// Write the RGB Data
	if (!WriteFile(hFile, pBitmapBits, bmpInfoHeader.biSizeImage, &dwBytesWritten, NULL))
	{
		CloseHandle(hFile);
		return E_FAIL;
	}

	// Close the file
	CloseHandle(hFile);
	return S_OK;
}