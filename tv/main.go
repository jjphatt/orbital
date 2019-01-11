package main

import (
  "net"
  "fmt"
  "image"
  "encoding/json"
  "strconv"
  "os"
  "github.com/fogleman/gg"
  "github.com/aarzilli/nucular"
)

// Here's the struct we unmarshal JSON data to.
// Note the field tags we use to translate lower-case JSON keys to 
// (Go exposed) upper-case fields.
type Datum struct {
  Name string `json:"name"`
  Time float64 `json:"time"`
  Data [16]float64 `json:"data"`
}

func CheckError(err error) {
  if err != nil {
    fmt.Println(err)
  }
}

// This goroutine receives UDP datagrams with JSON data and send them to the 
// given channel to be picked up by the Display goroutine.
func Receive(out chan<- Datum) {
  // Listen on a port.
  addr, err := net.ResolveUDPAddr("udp",":9999")
  CheckError(err)
  conn, err := net.ListenUDP("udp", addr)
  CheckError(err)

  buf := make([]byte, 2000)
  for {
    // Read a UDP packet.
    n, _, err := conn.ReadFromUDP(buf)
    if err != nil {
      fmt.Println(err)
      continue
    }

    // Decode the JSON contents and copy them into place.
    var data Datum
    err = json.Unmarshal(buf[:n], &data)
    if err != nil {
      fmt.Println(err)
      continue
    }
    out <- data
  }
}

type Canvas struct {
  Name string
  Frame, Index int
  w, h int
  image *image.RGBA
  gfx *gg.Context
  tMin, tMax, fMin, fMax float64
  t []float64
  f []float64
}

func NewCanvas(name string, 
               width, height int, 
               timeWindow, maxAmplitude float64,
               dataIndex int) *Canvas {
  canvas := Canvas{
    Name: name,
    Frame: 0,
    Index: dataIndex,
    image: nil,
    gfx: nil,
    w: width,
    h: height,
    tMin: -timeWindow,
    tMax: 0.0,
    fMin: -0.5 * maxAmplitude,
    fMax: 0.5 * maxAmplitude,
  }

  // Set up the graphics context
  canvas.image = image.NewRGBA(image.Rect(0, 0, width, height))
  canvas.gfx = gg.NewContextForRGBA(canvas.image)
  // Clear the graphics context.
  canvas.gfx.SetRGB(0, 0, 0)
  canvas.gfx.Clear()

  return &canvas
}

// Adds data to the given canvas, readjusting limits and discarding data
// that has fallen outside of the time window.
func (canvas *Canvas) AddData(data Datum) {
  // Rejigger the limits.
  tWin := canvas.tMax - canvas.tMin
  canvas.tMax = data.Time
  canvas.tMin = data.Time - tWin

  // Append new point.
  canvas.t = append(canvas.t, data.Time)
  canvas.f = append(canvas.f, data.Data[canvas.Index])

  // Remove old data outside of the window.
  for (canvas.t[0] < canvas.tMin) {
    canvas.t = canvas.t[1:]
    canvas.f = canvas.f[1:]
  }
}

func (canvas *Canvas) Draw(w *nucular.Window) {
  if len(canvas.t) > 0 {
    tWin := canvas.tMax - canvas.tMin
    x := ((canvas.t[0] - canvas.tMin) / tWin) * float64(canvas.w)
    y := float64(canvas.h)/2 + (canvas.f[0] / canvas.fMax) * float64(canvas.h/2)
    fmt.Printf("(x, y) = (%g, %g)\n", x, y);
    canvas.gfx.MoveTo(x, y)

    for i := 1; i < len(canvas.f); i++ {
      x = ((canvas.t[i] - canvas.tMin) / tWin) * float64(canvas.w)
      y = float64(canvas.h)/2 + (canvas.f[i] / canvas.fMax) * float64(canvas.h/2)
      canvas.gfx.LineTo(x, y)
    }

    w.Image(canvas.image)
  }
}

var win nucular.MasterWindow

// This goroutine pulls data out of the given channel and updates the master
// window.
func Display(canvas *Canvas, in <-chan Datum) {
  for {
    // Get the data from the channel.
    data := <-in

    // Add the new data to the canvas.
    canvas.AddData(data)
    fmt.Printf("Name: %s\nTime: %g\nData: (%g, %g, %g)\n", 
               data.Name, data.Time, data.Data[0], data.Data[1], data.Data[2])
    win.Changed();
  }
}

func main() {
  if len(os.Args) < 4 {
    fmt.Println("Usage: orbital.tv <timeWindow> <maxAmplitude> <dataIndex>")
    return
  }

  // Set up a new canvas.
  timeWindow, err := strconv.ParseFloat(os.Args[1], 64)
  CheckError(err)
  maxAmplitude, err := strconv.ParseFloat(os.Args[2], 64)
  CheckError(err)
  index, err := strconv.Atoi(os.Args[3])
  CheckError(err)
  canvas := NewCanvas("orbital", 1024, 768, timeWindow, maxAmplitude, index)

  // Set up the master window.
  win = nucular.NewMasterWindow(0, "orbital.tv", canvas.Draw)

  // Set up a channel and start some goroutines.
  ch := make(chan Datum)
  go Receive(ch)
  go Display(canvas, ch)

  win.Main()
}
