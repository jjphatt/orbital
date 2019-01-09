package main

import (
  "net"
  "fmt"
  "encoding/json"
  "bufio"
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
  gfx *gg.Context
  tMin, tMax, fMin, fMax float64
  data Datum
}

func NewCanvas(name string, 
               width, height int, 
               timeWindow, maxAmplitude float64,
               dataIndex int) *Canvas {
  canvas := Canvas{
    Name: name,
    Frame: 0,
    Index: dataIndex,
    gfx: gg.NewContext(width, height),
    tMin: -timeWindow,
    tMax: 0.0,
    fMin: -0.5 * maxAmplitude/2.0,
    fMax: 0.5 * maxAmplitude,
  }

  // Clear the graphics context.
  canvas.gfx.SetRGB(0, 0, 0)
  canvas.gfx.Clear()

  return &canvas
}

func (canvas *Canvas) Draw(w *nucular.Window) {
/*
  // Scale the points.
  x := data.Time /
  f := data.Data[canvas.Index]

  // If this is our first frame, Set our initial point.
  // Otherwise just connect the new dot.
  if canvas.Frame == 0 {
    canvas.gfx.MoveTo(x, y)
  } else {
    canvas.gfx.LineTo(x, y)
  }

  canvas.win.Image(canvas.gfx.Image())
  canvas.win.Changed()
*/
}

var win nucular.MasterWindow

func Display(canvas *Canvas, in <-chan Datum) {
  for {
    canvas.data = <-in
    fmt.Printf("Name: %s\nTime: %g\nData: (%g, %g, %g)\n", 
               canvas.data.Name, canvas.data.Time, canvas.data.Data[0], canvas.data.Data[1], canvas.data.Data[2])
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
  go win.Main()

  // Set up a channel and start some goroutines.
  ch := make(chan Datum)
  go Receive(ch)
  go Display(canvas, ch)

  // Break out of the loop when Escape is pressed.
  reader := bufio.NewReader(os.Stdin)
  for {
    char, _, _ := reader.ReadRune()
    if (char == 27) {
      break;
    }
  }
}
