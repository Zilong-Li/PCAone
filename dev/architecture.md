# Architecture

![PCAone class diagram](architecture.png)

The class diagram above is drawn with [PlantUML](https://plantuml.com) (1.2024.4)
from the source below. To redraw it after editing, run
`plantuml dev/architecture.md`: PlantUML picks up the `@startuml` block in this
file and writes `architecture.png` next to it.

```plantuml
@startuml architecture
skinparam dpi 300
skinparam defaultFontName FiraCode Nerd Font
skinparam groupInheritance 2
skinparam packageStyle rectangle
hide empty fields
hide empty members
abstract class Data {
  -Matrix A
  +void {abstract} read_all()
  +void {abstract} read_block()
}
Data <|-- FilePlink
Data <|-- FilePgen
Data <|-- FileBeagle
Data <|-- FileBgen
Data <|-- FileCSV
class FilePlink {
  +normalization()
}
class FilePgen {
  +normalization()
}
class FileBeagle {
  +normalization()
}
class FileBgen {
  +normalization()
}
class FileCSV {
  +normalization()
}
package "PCA methods" as PCA {
  class PCAone {
    +InCore
    +OutOfCore
  }
  class IRAM {
    +InCore
    +OutOfCore
  }
  class RSVD {
    +InCore
    +OutOfCore
  }
  class FullSVD {
    +InCore
  }
}

class "PCA Related" as PCAem {
  +LD
  +HWE
  +Projection
  +Selection
  +EMU
  +PCAngsd
}
FullSVD <-- Data: <:astonished:>
IRAM <-- Data: <:smile:>
RSVD <-- Data: <:nerd_face:>
PCAone <-- Data: <:sunglasses:> 
Data -left[dashed]-> PCAem: <:boom:> 
PCA -left-> PCAem
@enduml
```
