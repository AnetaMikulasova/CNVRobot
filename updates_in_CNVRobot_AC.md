Updates made in CNVRobot

1. Keeping unique controls in GATK CreateReadCountPanelOfNormals command

In cases where the same control samples appear multiple times in the control list, CreateReadCountPanelOfNormals fails because duplicate controls are included more than once. To prevent this, sort -u has been added to the command when creating a list of controls to retain only unique entries, eliminating duplicates and ensuring successful execution.